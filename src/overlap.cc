/*
    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "vdjsearch.h"

static struct db * d;
static unsigned int set1_longestsequence = 0;
static uint64_t set1_sequences = 0;
static uint64_t set1_residues = 0;
static uint64_t set1_samples = 0;
static uint64_t * set1_sample_size = nullptr;
static double * set1_sample_freq = nullptr;

static struct db * d2;
static unsigned int set2_longestsequence = 0;
static uint64_t set2_sequences = 0;
static uint64_t set2_residues = 0;
static uint64_t set2_samples = 0;
static uint64_t * set2_sample_size = nullptr;
static double * set2_sample_freq = nullptr;

static pthread_mutex_t network_mutex;
static uint64_t network_amp = 0;
static struct bloom_s * bloom_a = nullptr; // Bloom filter for sequences
static double * sample_matrix = nullptr;
static hashtable_s * hashtable = nullptr;

//#define AVOID_DUPLICATES 1

const uint64_t CHUNK = 1000;

void hash_insert(uint64_t amp)
{
  /* set 2 */
  /* find the first empty bucket */
  uint64_t hash = db_gethash(d2, amp);
  uint64_t j = hash_getindex(hashtable, hash);
  while (hash_is_occupied(hashtable, j))
    j = hash_getnextindex(hashtable, j);

  hash_set_occupied(hashtable, j);
  hash_set_value(hashtable, j, hash);
  hash_set_data(hashtable, j, amp);
  bloom_set(bloom_a, hash);
}

int set1_compare_by_sample_name(const void * a, const void * b)
{
  const unsigned int * x = (const unsigned int *) a;
  const unsigned int * y = (const unsigned int *) b;
  return strcmp(db_getsamplename(d, *x), db_getsamplename(d, *y));
}

int set2_compare_by_sample_name(const void * a, const void * b)
{
  const unsigned int * x = (const unsigned int *) a;
  const unsigned int * y = (const unsigned int *) b;
  return strcmp(db_getsamplename(d2, *x), db_getsamplename(d2, *y));
}

void find_variant_matches(uint64_t seed,
                          var_s * var,
                          double * sample_matrix,
                          struct bloom_s * bloom_d)
{
  /* seed in set 1, amp in set 2 */

  /* compute hash and corresponding hash table index */

  uint64_t j = hash_getindex(hashtable, var->hash);

  /* find matching buckets */

  while (hash_is_occupied(hashtable, j))
    {
      if (hash_compare_value(hashtable, j, var->hash))
        {
          uint64_t amp = hash_get_data(hashtable, j);

          /* double check that everything matches */

          unsigned int seed_v_gene = db_get_v_gene(d, seed);
          unsigned int seed_j_gene = db_get_j_gene(d, seed);

          unsigned int amp_v_gene = db_get_v_gene(d2, amp);
          unsigned int amp_j_gene = db_get_j_gene(d2, amp);

          if (opt_ignore_genes ||
              ((seed_v_gene == amp_v_gene) && (seed_j_gene == amp_j_gene)))
            {
              unsigned char * seed_sequence
                = (unsigned char *) db_getsequence(d, seed);
              unsigned int seed_seqlen
                = db_getsequencelen(d, seed);
              unsigned char * amp_sequence
                = (unsigned char *) db_getsequence(d2, amp);
              unsigned int amp_seqlen
                = db_getsequencelen(d2, amp);

              if (check_variant(seed_sequence, seed_seqlen,
                                var,
                                amp_sequence, amp_seqlen))
                {
                  unsigned int i = db_getsampleno(d, seed);
                  unsigned int j = db_getsampleno(d2, amp);
                  double f = db_get_freq(d, seed);
                  double g = db_get_freq(d2, amp);
                  sample_matrix[set2_samples * i + j] += f * g;

                  if (bloom_d)
                    bloom_set(bloom_d, var->hash);
                }
            }
        }
      j = hash_getnextindex(hashtable, j);
    }
}

void process_variants(uint64_t seed,
                      var_s * variant_list,
                      double * sample_matrix)
{
  unsigned int variant_count = 0;
  unsigned char * sequence = (unsigned char *) db_getsequence(d, seed);
  unsigned int seqlen = db_getsequencelen(d, seed);
  uint64_t hash = db_gethash(d, seed);
  uint64_t v_gene = db_get_v_gene(d, seed);
  uint64_t j_gene = db_get_j_gene(d, seed);

  generate_variants(hash,
                    sequence, seqlen, v_gene, j_gene,
                    variant_list, & variant_count);

  for(unsigned int i = 0; i < variant_count; i++)
    {
      var_s * var = variant_list + i;
      if (bloom_get(bloom_a, var->hash))
        {
          find_variant_matches(seed, var, sample_matrix, nullptr);
        }
    }
}

void process_variants_avoid_duplicates(uint64_t seed,
                                       var_s * variant_list,
                                       double * sample_matrix,
                                       struct bloom_s * bloom_d)
{
  unsigned int variant_count = 0;
  unsigned char * sequence = (unsigned char *) db_getsequence(d, seed);
  unsigned int seqlen = db_getsequencelen(d, seed);
  uint64_t hash = db_gethash(d, seed);
  uint64_t v_gene = db_get_v_gene(d, seed);
  uint64_t j_gene = db_get_j_gene(d, seed);

  generate_variants(hash,
                    sequence, seqlen, v_gene, j_gene,
                    variant_list, & variant_count);

  bloom_zap(bloom_d);

  for(unsigned int i = 0; i < variant_count; i++)
    {
      var_s * var = variant_list + i;

      if (bloom_get(bloom_a, var->hash))
        {

          /* Check for potential duplicate variants due to limitations */
          /* in the variant generation algorithm.                      */
          /* We only care about those duplicate that match the target. */

          bool dup = false;

          if (bloom_get(bloom_d, var->hash))
            {
              /* check if there is a real duplicate */

              printf("Potential duplicate variant!\n");

              unsigned char * seq1 = (unsigned char*) xmalloc(seqlen + 3);
              unsigned char * seq2 = (unsigned char*) xmalloc(seqlen + 3);

              for (unsigned int j = 0; j < i; j++)
                {
                  struct var_s * v = variant_list + j;

                  if (var->hash == v->hash)
                    {
                      printf("Likely duplicate variant!\n");
                      unsigned int seq1len, seq2len;
                      generate_variant_sequence(sequence,
                                                seqlen,
                                                var,
                                                seq1,
                                                & seq1len);

                      generate_variant_sequence(sequence,
                                                seqlen,
                                                var,
                                                seq2,
                                                & seq2len);

                      if ((seq1len == seq2len) &&
                          ! memcmp(seq1, seq2, seq1len))
                        {
                          /* we have a true duplicate variant */
                          printf("Real duplicate variant!\n");
                          dup = true;
                          break;
                        }
                    }
                }

              xfree(seq1);
              xfree(seq2);
            }

          if (!dup)
            find_variant_matches(seed, var, sample_matrix, bloom_d);
        }
    }
}

void sim_thread(int64_t t)
{
  (void) t;

  uint64_t maxvar = max_variants(set1_longestsequence);

  struct var_s * variant_list = static_cast<struct var_s *>
    (xmalloc(maxvar * sizeof(struct var_s)));

  double * sample_matrix_local = static_cast<double *>
    (xmalloc(set1_samples * set2_samples * sizeof(double)));

  for(uint64_t k = 0; k < set1_samples * set2_samples; k++)
    sample_matrix_local[k] = 0;

#ifdef AVOID_DUPLICATES
  /* init bloom filter for duplicates */
  uint64_t bloomsize = 1;
  while (bloomsize < maxvar)
    bloomsize *= 2;
  struct bloom_s * bloom_d = bloom_init(bloomsize);
#endif

  pthread_mutex_lock(&network_mutex);

  while (network_amp < set1_sequences)
    {
      uint64_t firstseed = network_amp;
      network_amp += CHUNK;
      if (network_amp > set1_sequences)
        network_amp = set1_sequences;
      progress_update(network_amp);
      pthread_mutex_unlock(&network_mutex);

      uint64_t chunksize = network_amp - firstseed;

      /* process chunksize sequences starting at seed */
      for (uint64_t z = 0; z < chunksize; z++)
        {
          uint64_t seed = firstseed + z;

#ifdef AVOID_DUPLICATES
          process_variants_avoid_duplicates(seed, variant_list, sample_matrix_local, bloom_d);
#else
          process_variants(seed, variant_list, sample_matrix_local);
#endif
        }

      /* lock mutex and update global sample_matrix */
      pthread_mutex_lock(&network_mutex);
    }

  /* update global sample_matrix */
  for(uint64_t k = 0; k < set1_samples * set2_samples; k++)
    sample_matrix[k] += sample_matrix_local[k];

  pthread_mutex_unlock(&network_mutex);

#ifdef AVOID_DUPLICATES
  bloom_exit(bloom_d);
#endif

  xfree(sample_matrix_local);
  xfree(variant_list);
}

void overlap(char * set1_filename, char * set2_filename)
{
  /* find overlaps between repertoires of samples */

  db_init();


  /**** Set 1 ****/

  fprintf(logfile, "Immune receptor repertoire set 1\n");

  d = db_create();
  db_read(d, set1_filename);

  set1_longestsequence = db_getlongestsequence(d);
  set1_sequences = db_getsequencecount(d);
  set1_samples = db_getsamplecount(d);
  set1_residues = db_getresiduescount(d);

  fprintf(logfile, "\n");

  /* determine number of sequences in each of the samples (Set 1) */

  set1_sample_size = static_cast<uint64_t *>
    (xmalloc(sizeof(uint64_t) * set1_samples));
  for (unsigned int s = 0; s < set1_samples ; s++)
    set1_sample_size[s] = 0;
  set1_sample_freq = static_cast<double *>
    (xmalloc(sizeof(double) * set1_samples));
  for (unsigned int s = 0; s < set1_samples ; s++)
    set1_sample_freq[s] = 0.0;
  for (uint64_t i = 0; i < set1_sequences ; i++)
    {
      unsigned int s = db_getsampleno(d, i);
      set1_sample_size[s]++;
      set1_sample_freq[s] += db_get_freq(d, i);
    }

  /* set 1 : sort samples alphanumerically for display */

  unsigned int * set1_lookup_sample =
    (unsigned int *) xmalloc(sizeof(unsigned int) * set1_samples);
  for (unsigned int i = 0; i < set1_samples; i++)
    set1_lookup_sample[i] = i;
  qsort(set1_lookup_sample,
        set1_samples,
        sizeof(unsigned int),
        set1_compare_by_sample_name);

  /* list of samples in set 1 */

  fprintf(logfile, "#no\tseqs\tfreq\tsample\n");
  uint64_t sum_size = 0;
  double sum_freq = 0.0;
  for (unsigned int i = 0; i < set1_samples; i++)
    {
      unsigned int s = set1_lookup_sample[i];
      if (opt_ignore_frequency)
        fprintf(logfile, "%u\t%7" PRIu64 "\t%7.0lf\t%s\n",
                i+1,
                set1_sample_size[s],
                set1_sample_freq[s],
                db_getsamplename(d, s));
      else
        fprintf(logfile, "%u\t%7" PRIu64 "\t%7.5lf\t%s\n",
                i+1,
                set1_sample_size[s],
                set1_sample_freq[s],
                db_getsamplename(d, s));
      sum_size += set1_sample_size[s];
      sum_freq += set1_sample_freq[s];
    }
  if (opt_ignore_frequency)
    fprintf(logfile, "Sum\t%" PRIu64 "\t%7.0lf\n", sum_size, sum_freq);
  else
    fprintf(logfile, "Sum\t%" PRIu64 "\t%7.5lf\n", sum_size, sum_freq);
  fprintf(logfile, "\n");


  /**** Set 2 ****/

  fprintf(logfile, "Immune receptor repertoire set 2\n");

  d2 = db_create();
  db_read(d2, set2_filename);

  set2_longestsequence = db_getlongestsequence(d2);
  set2_sequences = db_getsequencecount(d2);
  set2_samples = db_getsamplecount(d2);
  set2_residues = db_getresiduescount(d2);

  fprintf(logfile, "\n");

  /* determine number of sequences in each of the samples (Set 2) */

  set2_sample_size = static_cast<uint64_t *>
    (xmalloc(sizeof(uint64_t) * set2_samples));
  for (unsigned int t = 0; t < set2_samples ; t++)
    set2_sample_size[t] = 0;
  set2_sample_freq = static_cast<double *>
    (xmalloc(sizeof(double) * set2_samples));
  for (unsigned int t = 0; t < set2_samples ; t++)
    set2_sample_freq[t] = 0.0;
  for (uint64_t j = 0; j < set2_sequences ; j++)
    {
      unsigned int t = db_getsampleno(d2, j);
      set2_sample_size[t]++;
      set2_sample_freq[t] += db_get_freq(d2, j);
    }

  /* set 2 : sort samples alphanumerically for display */

  unsigned int * set2_lookup_sample =
    (unsigned int *) xmalloc(sizeof(unsigned int) * set2_samples);
  for (unsigned int j = 0; j < set2_samples; j++)
    set2_lookup_sample[j] = j;
  qsort(set2_lookup_sample,
        set2_samples,
        sizeof(unsigned int),
        set2_compare_by_sample_name);

  /* list of samples in set 2 */

  sum_size = 0;
  sum_freq = 0.0;
  fprintf(logfile, "#no\tseqs\tfreq\tsample\n");
  for (unsigned int j = 0; j < set2_samples; j++)
    {
      unsigned int t = set2_lookup_sample[j];
      if (opt_ignore_frequency)
        fprintf(logfile, "%u\t%7" PRIu64 "\t%7.0lf\t%s\n",
                j+1,
                set2_sample_size[t],
                set2_sample_freq[t],
                db_getsamplename(d2, t));
      else
        fprintf(logfile, "%u\t%7" PRIu64 "\t%7.5lf\t%s\n",
                j+1,
                set2_sample_size[t],
                set2_sample_freq[t],
                db_getsamplename(d2, t));
      sum_size += set2_sample_size[t];
      sum_freq += set2_sample_freq[t];
    }
  if (opt_ignore_frequency)
    fprintf(logfile, "Sum\t%" PRIu64 "\t%7.0lf\n", sum_size, sum_freq);
  else
    fprintf(logfile, "Sum\t%" PRIu64 "\t%7.5lf\n", sum_size, sum_freq);
  fprintf(logfile, "\n");

  unsigned int overall_longest = MAX(set1_longestsequence,
                                     set2_longestsequence);

  zobrist_init(overall_longest + 2,
               db_get_v_gene_count(),
               db_get_j_gene_count());

  fprintf(logfile, "Unique V genes:    %" PRIu64 "\n",
          db_get_v_gene_count());

  fprintf(logfile, "Unique J genes:    %" PRIu64 "\n",
          db_get_j_gene_count());

  /* compute hashes for each sequence in database */

  db_hash(d);
  db_hash(d2);

  /* compute hash for all sequences and store them in a hash table */
  /* use an additional bloom filter for increased speed */
  /* hashing into hash table & bloom filter */

  hashtable = hash_init(set2_sequences);
  bloom_a = bloom_init(hash_get_tablesize(hashtable) * 2);
  progress_init("Hashing sequences:", set2_sequences);
  for(uint64_t i=0; i < set2_sequences; i++)
    {
      hash_insert(i);
      progress_update(i);
    }
  progress_done();

  /* allocate matrix of sample1 x sample2 counts */

  sample_matrix = static_cast<double *>
    (xmalloc(sizeof(double) * set1_samples * set2_samples));

  for(unsigned int s = 0; s < set1_samples; s++)
    for(unsigned int t = 0; t < set2_samples; t++)
      sample_matrix[set2_samples * s + t] = 0;

  /* compare all sequences */

  pthread_mutex_init(&network_mutex, nullptr);
  progress_init("Analysing:        ", set1_sequences);

  if (opt_threads == 1)
    {
      sim_thread(0);
    }
  else
    {
      ThreadRunner * sim_tr = new ThreadRunner(static_cast<int>(opt_threads),
                                               sim_thread);
      sim_tr->run();
      delete sim_tr;
    }

  progress_done();
  pthread_mutex_destroy(&network_mutex);

  /* dump similarity matrix */

  unsigned int x = 0;
  progress_init("Writing results:  ", set1_sequences * set2_sequences);
  if (opt_alternative)
    {
      for (unsigned int i = 0; i < set1_samples; i++)
        {
          unsigned int s = set1_lookup_sample[i];
          for (unsigned int j = 0; j < set2_samples; j++)
            {
              unsigned int t = set2_lookup_sample[j];
              fprintf(outfile,
                      "%s\t%s\t%15.9le\n",
                      db_getsamplename(d, s),
                      db_getsamplename(d2, t),
                      sample_matrix[set2_samples * s + t]);
              x++;
            }
        }
    }
  else
    {
      fprintf(outfile, "#sample");
      for (unsigned int j = 0; j < set2_samples; j++)
        fprintf(outfile, "\t%s", db_getsamplename(d2, set2_lookup_sample[j]));
      fprintf(outfile, "\n");
      for (unsigned int i = 0; i < set1_samples; i++)
        {
          unsigned int s = set1_lookup_sample[i];
          fprintf(outfile, "%s", db_getsamplename(d, s));
          for (unsigned int j = 0; j < set2_samples; j++)
            {
              unsigned int t = set2_lookup_sample[j];
              fprintf(outfile, "\t%15.9le", sample_matrix[set2_samples * s + t]);
              x++;
              progress_update(x);
            }
          fprintf(outfile, "\n");
        }
    }
  progress_done();

  fprintf(logfile, "\n");

  if (sample_matrix)
    xfree(sample_matrix);

  bloom_exit(bloom_a);

  hash_exit(hashtable);

  db_free(d2);
  db_free(d);
  db_exit();

  zobrist_exit();

  xfree(set2_lookup_sample);
  xfree(set2_sample_freq);
  xfree(set2_sample_size);

  xfree(set1_lookup_sample);
  xfree(set1_sample_freq);
  xfree(set1_sample_size);
}
