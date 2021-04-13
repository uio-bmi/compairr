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
static unsigned int set1_sequences = 0;
static uint64_t set1_residues = 0;
static uint64_t set1_samples = 0;
static uint64_t * set1_sample_size = nullptr;
static double * set1_sample_freq = nullptr;

static struct db * d2;
static unsigned int set2_longestsequence = 0;
static unsigned int set2_sequences = 0;
static uint64_t set2_residues = 0;
static uint64_t set2_samples = 0;
static uint64_t * set2_sample_size = nullptr;
static double * set2_sample_freq = nullptr;

static pthread_mutex_t network_mutex;
static unsigned int network_amp = 0;
static struct bloom_s * bloom_a = nullptr; // Bloom filter for sequences
static double * sample_matrix = nullptr;
static hashtable_s * hashtable = nullptr;

//#define COUNTMATCHES
#ifdef COUNTMATCHES
static uint64_t matches = 0;
#endif

//#define AVOID_DUPLICATES 1

const unsigned int CHUNK = 1000;

inline void hash_insert(unsigned int amp)
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

void find_variant_matches(unsigned int seed,
                          var_s * var,
                          double * sample_hits,
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
          unsigned int amp = hash_get_data(hashtable, j);

          /* double check that everything matches */

          unsigned int seed_v_gene = db_get_v_gene(d, seed);
          unsigned int seed_d_gene = db_get_d_gene(d, seed);

          unsigned int amp_v_gene = db_get_v_gene(d2, amp);
          unsigned int amp_d_gene = db_get_d_gene(d2, amp);

          if (opt_ignore_genes ||
              ((seed_v_gene == amp_v_gene) && (seed_d_gene == amp_d_gene)))
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
                  double g = db_get_freq(d2, amp);
                  sample_hits[db_getsampleno(d2, amp)] += g;

#ifdef COUNTMATCHES
                  matches++;
#endif

                  if (bloom_d)
                    bloom_set(bloom_d, var->hash);
                }
            }
        }
      j = hash_getnextindex(hashtable, j);
    }
}

void process_variants(unsigned int seed,
                      var_s * variant_list,
                      double * sample_hits)
{
  unsigned int variant_count = 0;
  unsigned char * sequence = (unsigned char *) db_getsequence(d, seed);
  unsigned int seqlen = db_getsequencelen(d, seed);
  uint64_t hash = db_gethash(d, seed);
  uint64_t v_gene = db_get_v_gene(d, seed);
  uint64_t d_gene = db_get_d_gene(d, seed);

  generate_variants(hash,
                    sequence, seqlen, v_gene, d_gene,
                    variant_list, & variant_count);

  for(unsigned int i = 0; i < variant_count; i++)
    {
      var_s * var = variant_list + i;
      if (bloom_get(bloom_a, var->hash))
        {
          find_variant_matches(seed, var, sample_hits, nullptr);
        }
    }
}

void process_variants_avoid_duplicates(unsigned int seed,
                                       var_s * variant_list,
                                       double * sample_hits,
                                       struct bloom_s * bloom_d)
{
  unsigned int variant_count = 0;
  unsigned char * sequence = (unsigned char *) db_getsequence(d, seed);
  unsigned int seqlen = db_getsequencelen(d, seed);
  uint64_t hash = db_gethash(d, seed);
  uint64_t v_gene = db_get_v_gene(d, seed);
  uint64_t d_gene = db_get_d_gene(d, seed);

  generate_variants(hash,
                    sequence, seqlen, v_gene, d_gene,
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
            find_variant_matches(seed, var, sample_hits, bloom_d);
        }
    }
}

void sim_thread(int64_t t)
{
  (void) t;

  uint64_t maxvar = max_variants(set1_longestsequence);

  struct var_s * variant_list = static_cast<struct var_s *>
    (xmalloc(maxvar * sizeof(struct var_s)));

  double * sample_hits = static_cast<double *>
    (xmalloc(set2_samples * sizeof(double)));

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
      unsigned int firstseed = network_amp;
      network_amp += CHUNK;
      if (network_amp > set1_sequences)
        network_amp = set1_sequences;
      progress_update(network_amp);
      pthread_mutex_unlock(&network_mutex);

      unsigned int chunksize = network_amp - firstseed;

      /* process chunksize sequences starting at seed */
      for (unsigned int z = 0; z < chunksize; z++)
        {
          unsigned int seed = firstseed + z;
          uint64_t i = db_getsampleno(d, seed);
          double f = db_get_freq(d, seed);

          for(uint64_t j = 0; j < set2_samples; j++)
            sample_hits[j] = 0;

#ifdef AVOID_DUPLICATES
          process_variants_avoid_duplicates(seed, variant_list, sample_hits, bloom_d);
#else
          process_variants(seed, variant_list, sample_hits);
#endif
          uint64_t base = set2_samples * i;
          for(uint64_t j = 0; j < set2_samples; j++)
            sample_matrix_local[base + j] += sample_hits[j] * f;
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
  xfree(sample_hits);
  xfree(variant_list);
}

void overlap(char * set1_filename, char * set2_filename)
{
  /* find overlaps between repertoires of samples */

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
  for (unsigned int i = 0; i < set1_sequences ; i++)
    {
      unsigned int s = db_getsampleno(d, i);
      set1_sample_size[s]++;
      set1_sample_freq[s] += db_get_freq(d, i);
    }

  /* set 1 : sort samples alphanumerically for display */

  unsigned int * set1_lookup_sample =
    (unsigned int *) xmalloc(sizeof(uint64_t) * set1_samples);
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
  for (unsigned int j = 0; j < set2_sequences ; j++)
    {
      unsigned int t = db_getsampleno(d2, j);
      set2_sample_size[t]++;
      set2_sample_freq[t] += db_get_freq(d2, j);
    }

  /* set 2 : sort samples alphanumerically for display */

  unsigned int * set2_lookup_sample =
    (unsigned int *) xmalloc(sizeof(uint64_t) * set2_samples);
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
               db_get_d_gene_count());

  fprintf(logfile, "Unique v_genes:    %" PRIu64 "\n",
          db_get_v_gene_count());

  fprintf(logfile, "Unique d_genes:    %" PRIu64 "\n",
          db_get_d_gene_count());

  /* compute hashes for each sequence in database */

  db_hash(d);
  db_hash(d2);

  /* compute hash for all sequences and store them in a hash table */
  /* use an additional bloom filter for increased speed */
  /* hashing into hash table & bloom filter */

  hashtable = hash_init(set2_sequences);
  bloom_a = bloom_init(hash_get_tablesize(hashtable) * 2);
  progress_init("Hashing sequences:", set2_sequences);
  for(unsigned int i=0; i < set2_sequences; i++)
    {
      hash_insert(i);
      progress_update(i);
    }
  progress_done();

  /* allocate matrix of sample x sample counts */

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
  pthread_mutex_destroy(&network_mutex);
  progress_done();

  fprintf(logfile, "\n");

  /* dump similarity matrix */

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
          fprintf(outfile, "\t%7.1le", sample_matrix[set2_samples * s + t]);
        }
      fprintf(outfile, "\n");
    }

  bloom_exit(bloom_a);

  xfree(sample_matrix);

  xfree(set1_sample_size);
  xfree(set2_sample_size);

  xfree(set1_lookup_sample);
  xfree(set2_lookup_sample);

  zobrist_exit();

  db_free(d);
  db_free(d2);
}
