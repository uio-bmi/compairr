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

#include "compairr.h"

static struct db * d1;
static unsigned int set1_longestsequence = 0;
static uint64_t set1_sequences = 0;
static uint64_t set1_residues = 0;
static uint64_t set1_repertoires = 0;
static uint64_t * set1_repertoire_size = nullptr;
static uint64_t * set1_repertoire_count = nullptr;
static double * set1_repertoire_sq_count = nullptr;
static unsigned int * set1_lookup_repertoire = nullptr;

static struct db * d2;
static unsigned int set2_longestsequence = 0;
static uint64_t set2_sequences = 0;
static uint64_t set2_residues = 0;
static uint64_t set2_repertoires = 0;
static uint64_t * set2_repertoire_size = nullptr;
static uint64_t * set2_repertoire_count = nullptr;
static double * set2_repertoire_sq_count = nullptr;
static unsigned int * set2_lookup_repertoire = nullptr;

typedef double m_val_t;

static pthread_mutex_t pairs_mutex;
static pthread_mutex_t network_mutex;
static uint64_t network_progress = 0;
static struct bloom_s * bloom_a = nullptr; // Bloom filter for sequences
static m_val_t * repertoire_matrix = nullptr;
static hashtable_s * hashtable = nullptr;

static uint64_t all_matches = 0;

typedef struct pair_s
{
  uint64_t seq[2];
} pair_t;

const uint64_t CHUNK = 1000;
const char * empty_string = "";

static inline bool hash_insert(struct db * d,
                               hashtable_s * ht,
                               struct bloom_s * b,
                               uint64_t seed)
{
  bool dup = false;

  /* find the first empty bucket */
  uint64_t hash = db_gethash(d, seed);
  uint64_t j = hash_getindex(ht, hash);
  while (hash_is_occupied(ht, j))
    {
#if 1
      if (((b == nullptr) || bloom_get(b, hash)) &&
          (hash_compare_value(ht, j, hash)))
#else
      if (hash_compare_value(ht, j, hash))
#endif
        {
          uint64_t hit = hash_get_data(ht, j);

          /* check repertoire id match */
          unsigned int seed_rep_id = db_get_repertoire_id_no(d, seed);
          unsigned int hit_rep_id = db_get_repertoire_id_no(d, hit);

          if (seed_rep_id == hit_rep_id)
            {
              /* double check that everything matches */
              unsigned int seed_v_gene = db_get_v_gene(d, seed);
              unsigned int seed_j_gene = db_get_j_gene(d, seed);

              unsigned int hit_v_gene = db_get_v_gene(d, hit);
              unsigned int hit_j_gene = db_get_j_gene(d, hit);

              if (opt_ignore_genes ||
                  ((seed_v_gene == hit_v_gene) && (seed_j_gene == hit_j_gene)))
                {
                  unsigned char * seed_sequence
                    = (unsigned char *) db_getsequence(d, seed);
                  unsigned int seed_seqlen
                    = db_getsequencelen(d, seed);
                  unsigned char * hit_sequence
                    = (unsigned char *) db_getsequence(d, hit);
                  unsigned int hit_seqlen
                    = db_getsequencelen(d, hit);

                  if ((seed_seqlen == hit_seqlen) &&
                      ! memcmp(seed_sequence, hit_sequence, seed_seqlen))
                    {
                      dup = true;
                    }
                }
            }
        }
      j = hash_getnextindex(ht, j);
    }

  hash_set_occupied(ht, j);
  hash_set_value(ht, j, hash);
  hash_set_data(ht, j, seed);

  if (b)
    bloom_set(b, hash);

  return dup;
}

static int set1_compare_by_repertoire_name(const void * a, const void * b)
{
  const unsigned int * x = (const unsigned int *) a;
  const unsigned int * y = (const unsigned int *) b;
  return strcmp(db_get_repertoire_id(d1, *x), db_get_repertoire_id(d1, *y));
}

static int set2_compare_by_repertoire_name(const void * a, const void * b)
{
  const unsigned int * x = (const unsigned int *) a;
  const unsigned int * y = (const unsigned int *) b;
  return strcmp(db_get_repertoire_id(d2, *x), db_get_repertoire_id(d2, *y));
}

static inline m_val_t compute_score(uint64_t a, uint64_t b)
{
  if (opt_ignore_counts)
    return 1;
  else
    switch(opt_score_int)
      {
      case score_mh:
      case score_product:
        return (m_val_t)(a) * (m_val_t)(b);
      case score_ratio:
        return (m_val_t)(a) / (m_val_t)(b);
      case score_jaccard:
      case score_min:
        return MIN(a, b);
      case score_max:
        return MAX(a, b);
      case score_mean:
        return ((m_val_t)(a) + (m_val_t)(b)) / 2;
      default:
        fatal("Internal error");
      }
}

static void find_variant_matches(uint64_t seed,
                                 var_s * var,
                                 m_val_t * repertoire_matrix,
                                 uint64_t * pairs_alloc,
                                 uint64_t * pairs_count,
                                 struct pair_s * * pairs_list)
{
  /* compute hash and corresponding hash table index */

  uint64_t j = hash_getindex(hashtable, var->hash);

  /* find matching buckets */

  while (hash_is_occupied(hashtable, j))
    {
      if (hash_compare_value(hashtable, j, var->hash))
        {
          uint64_t hit = hash_get_data(hashtable, j);

          /* double check that everything matches */

          unsigned int seed_v_gene = db_get_v_gene(d1, seed);
          unsigned int seed_j_gene = db_get_j_gene(d1, seed);

          unsigned int hit_v_gene = db_get_v_gene(d2, hit);
          unsigned int hit_j_gene = db_get_j_gene(d2, hit);

          if (opt_ignore_genes ||
              ((seed_v_gene == hit_v_gene) && (seed_j_gene == hit_j_gene)))
            {
              unsigned char * seed_sequence
                = (unsigned char *) db_getsequence(d1, seed);
              unsigned int seed_seqlen
                = db_getsequencelen(d1, seed);
              unsigned char * hit_sequence
                = (unsigned char *) db_getsequence(d2, hit);
              unsigned int hit_seqlen
                = db_getsequencelen(d2, hit);

              if (check_variant(seed_sequence, seed_seqlen,
                                var,
                                hit_sequence, hit_seqlen))
                {
                  unsigned int i = db_get_repertoire_id_no(d1, seed);
                  unsigned int j = db_get_repertoire_id_no(d2, hit);
                  uint64_t f = db_get_count(d1, seed);
                  uint64_t g = db_get_count(d2, hit);

                  m_val_t s = compute_score(f, g);

		  if (! opt_no_matrix)
		    {
		      if (opt_matrix)
			{
			  repertoire_matrix[set2_repertoires * i + j] += s;
			}
		      else
			{
			  repertoire_matrix[set2_repertoires * seed + j] += s;
			}
		    }

                  all_matches++;

                  if (opt_pairs)
                    {
                      /* allocate more memory if needed */
                      if (*pairs_count >= *pairs_alloc)
                        {
                          * pairs_alloc = 2 * (* pairs_alloc);
                          * pairs_list = static_cast<struct pair_s *>
                            (xrealloc(* pairs_list,
                                      (*pairs_alloc) * sizeof(struct pair_s)));
                        }

                      struct pair_s p = { seed, hit };
                      (*pairs_list)[(*pairs_count)++] = p;
                    }
                }
            }
        }
      j = hash_getnextindex(hashtable, j);
    }
}

static void process_variants(uint64_t seed,
                             var_s * variant_list,
                             m_val_t * repertoire_matrix,
                             uint64_t * pairs_alloc,
                             uint64_t * pairs_count,
                             struct pair_s * * pairs_list)
{
  unsigned int variant_count = 0;
  unsigned char * sequence = (unsigned char *) db_getsequence(d1, seed);
  unsigned int seqlen = db_getsequencelen(d1, seed);
  uint64_t hash = db_gethash(d1, seed);
  uint64_t v_gene = db_get_v_gene(d1, seed);
  uint64_t j_gene = db_get_j_gene(d1, seed);

  generate_variants(hash,
                    sequence, seqlen, v_gene, j_gene,
                    variant_list, & variant_count);

  for(unsigned int i = 0; i < variant_count; i++)
    {
      var_s * var = variant_list + i;
      if (bloom_get(bloom_a, var->hash))
        {
          find_variant_matches(seed,
                               var,
                               repertoire_matrix,
                               pairs_alloc,
                               pairs_count,
                               pairs_list);
        }
    }
}

static void process_trad(uint64_t seed,
                         m_val_t * repertoire_matrix,
                         uint64_t * pairs_alloc,
                         uint64_t * pairs_count,
                         struct pair_s * * pairs_list)
{
  /* Only to be used with no indels (and d >= 3) */

  for (uint64_t hit = 0; hit < set2_sequences; hit++)
    {
      /* check if everything matches */

      unsigned int seed_v_gene = db_get_v_gene(d1, seed);
      unsigned int seed_j_gene = db_get_j_gene(d1, seed);

      unsigned int hit_v_gene = db_get_v_gene(d2, hit);
      unsigned int hit_j_gene = db_get_j_gene(d2, hit);

      if (opt_ignore_genes ||
          ((seed_v_gene == hit_v_gene) && (seed_j_gene == hit_j_gene)))
        {
          unsigned int seed_seqlen = db_getsequencelen(d1, seed);
          unsigned int hit_seqlen = db_getsequencelen(d2, hit);

          if (seed_seqlen == hit_seqlen)
            {
              unsigned char * seed_sequence
                = (unsigned char *) db_getsequence(d1, seed);
              unsigned char * hit_sequence
                = (unsigned char *) db_getsequence(d2, hit);

              if (seq_diff(seed_sequence, hit_sequence, seed_seqlen)
                  <= opt_differences)
              {
                unsigned int i = db_get_repertoire_id_no(d1, seed);
                unsigned int j = db_get_repertoire_id_no(d2, hit);
                uint64_t f = db_get_count(d1, seed);
                uint64_t g = db_get_count(d2, hit);

                m_val_t s = compute_score(f, g);

                if (! opt_no_matrix)
                  {
                    if (opt_matrix)
                      {
                        repertoire_matrix[set2_repertoires * i + j] += s;
                      }
                    else
                      {
                        repertoire_matrix[set2_repertoires * seed + j] += s;
                      }
                  }

                all_matches++;

                if (opt_pairs)
                  {
                    /* allocate more memory if needed */
                    if (*pairs_count >= *pairs_alloc)
                      {
                        * pairs_alloc = 2 * (* pairs_alloc);
                        * pairs_list = static_cast<struct pair_s *>
                          (xrealloc(* pairs_list,
                                    (*pairs_alloc) * sizeof(struct pair_s)));
                      }

                    struct pair_s p = { seed, hit };
                    (*pairs_list)[(*pairs_count)++] = p;
                  }
              }
            }
        }
    }
}

static void process_seq(uint64_t seed,
                        var_s * variant_list,
                        m_val_t * repertoire_matrix,
                        uint64_t * pairs_alloc,
                        uint64_t * pairs_count,
                        struct pair_s * * pairs_list)
{
  if (opt_differences > MAXDIFF_HASH)
    process_trad(seed, repertoire_matrix,
                 pairs_alloc, pairs_count, pairs_list);
  else
    process_variants(seed, variant_list, repertoire_matrix,
                     pairs_alloc, pairs_count, pairs_list);
}

static void sim_thread(int64_t t)
{
  (void) t;

  uint64_t maxvar = max_variants(set1_longestsequence);

  uint64_t pairs_alloc = 4 * CHUNK;
  uint64_t pairs_count = 0;

  struct pair_s * pairs_list = nullptr;
  if (opt_pairs)
    pairs_list = static_cast<struct pair_s *>
      (xmalloc(pairs_alloc * sizeof(struct pair_s)));

  struct var_s * variant_list = static_cast<struct var_s *>
    (xmalloc(maxvar * sizeof(struct var_s)));

  m_val_t * repertoire_matrix_local = nullptr;
  if (opt_threads > 1)
    {
      /* if multiple threads, create local matrix */

      if (! opt_no_matrix)
        {
          if (opt_matrix)
            {
              repertoire_matrix_local = static_cast<m_val_t *>
                (xmalloc(set1_repertoires * set2_repertoires * sizeof(m_val_t)));

              for(uint64_t k = 0; k < set1_repertoires * set2_repertoires; k++)
                repertoire_matrix_local[k] = 0;
            }
          else
            {
              repertoire_matrix_local = static_cast<m_val_t *>
                (xmalloc(set1_sequences * set2_repertoires * sizeof(m_val_t)));

              for(uint64_t k = 0; k < set1_sequences * set2_repertoires; k++)
                repertoire_matrix_local[k] = 0;
            }
        }

      pthread_mutex_lock(&network_mutex);
    }

  while (network_progress < set1_sequences)
    {
      uint64_t firstseed = network_progress;
      network_progress += CHUNK;
      if (network_progress > set1_sequences)
        network_progress = set1_sequences;
      progress_update(network_progress);
      uint64_t chunksize = network_progress - firstseed;

      if (opt_threads > 1)
        {
          pthread_mutex_unlock(&network_mutex);
        }

      /* process chunksize sequences starting at seed */

      for (uint64_t z = 0; z < chunksize; z++)
        {
          uint64_t seed = firstseed + z;
          process_seq(seed,
                      variant_list,
                      (opt_threads > 1 ?
                       repertoire_matrix_local :
                       repertoire_matrix),
                      & pairs_alloc,
                      & pairs_count,
                      & pairs_list);
        }

      if (opt_threads > 1)
        {
          pthread_mutex_lock(&network_mutex);
        }

      if (opt_pairs)
        {
          for (uint64_t i = 0; i < pairs_count; i++)
            {
              uint64_t a = pairs_list[i].seq[0];
              uint64_t b = pairs_list[i].seq[1];

              int rep_id_no1 = db_get_repertoire_id_no(d1, a);
              const char * rep_id1 = db_get_repertoire_id(d1, rep_id_no1);
              int64_t len1 = db_getsequencelen(d1, a);

              int rep_id_no2 = db_get_repertoire_id_no(d2, b);
              const char * rep_id2 = db_get_repertoire_id(d2, rep_id_no2);
              int64_t len2 = db_getsequencelen(d2, b);

              fprintf(pairsfile,
                      "%s\t%s\t%" PRIu64 "\t%s\t%s\t",
                      rep_id1,
                      db_get_sequence_id(d1, a),
                      db_get_count(d1, a),
                      db_get_v_gene_name(d1, a),
                      db_get_j_gene_name(d1, a));
              db_fprint_sequence(pairsfile, d1, a);
              if (opt_keep_columns)
                fprintf(pairsfile, "\t%s", db_get_keep_columns(d1, a));

              fprintf(pairsfile,
                      "\t%s\t%s\t%" PRIu64 "\t%s\t%s\t",
                      rep_id2,
                      db_get_sequence_id(d2, b),
                      db_get_count(d2, b),
                      db_get_v_gene_name(d2, b),
                      db_get_j_gene_name(d2, b));
              db_fprint_sequence(pairsfile, d2, b);
              if (opt_keep_columns)
                fprintf(pairsfile, "\t%s", db_get_keep_columns(d2, b));

              if (opt_distance)
                {
                  /* Compute Hamming distance if len1 = len2,
                     otherwise use Levenshtein distance of 1 (one indel) */
                  int64_t dist = 1;
                  if (len1 == len2)
                    dist = seq_diff((unsigned char *)db_getsequence(d1, a),
                                    (unsigned char *)db_getsequence(d2, b),
                                    len1);
                  fprintf(pairsfile, "\t%" PRId64, dist);
                }

              fprintf(pairsfile, "\n");
            }
          pairs_count = 0;
        }
    }

  if (opt_threads > 1)
    {
      /* update global repertoire_matrix */
      if (! opt_no_matrix)
        {
          if (opt_matrix)
            {
              for(uint64_t k = 0; k < set1_repertoires * set2_repertoires; k++)
                repertoire_matrix[k] += repertoire_matrix_local[k];
            }
          else
            {
              for(uint64_t k = 0; k < set1_sequences * set2_repertoires; k++)
                repertoire_matrix[k] += repertoire_matrix_local[k];
            }
        }

      pthread_mutex_unlock(&network_mutex);

      if (repertoire_matrix_local)
        xfree(repertoire_matrix_local);
      repertoire_matrix_local = nullptr;
    }

  xfree(variant_list);

  if (opt_pairs)
    xfree(pairs_list);
}

static void show_matrix_value(unsigned int s, unsigned int t)
{
  double SP, LX, LY, XY, MH;
  double SM, SA, SB, JI;
  double X;

  switch (opt_score_int)
    {
    case score_mh:
      /* Morisita-Horn index */
      /* Uses sum of product */
      SP = repertoire_matrix[set2_repertoires * s + t];
      LX = set1_repertoire_sq_count[s] /
        set1_repertoire_count[s] / set1_repertoire_count[s];
      LY = set2_repertoire_sq_count[t] /
        set2_repertoire_count[t] / set2_repertoire_count[t];
      XY = 1.0 * set1_repertoire_count[s]
        * set2_repertoire_count[t];
      MH = (2.0 * SP) / ((LX + LY) * XY);
      fprintf(outfile, "\t%.10lg", MH);
      break;

    case score_jaccard:
      /* Jaccard index */
      /* Uses sum of min */
      SM = repertoire_matrix[set2_repertoires * s + t];
      SA = set1_repertoire_count[s];
      SB = set2_repertoire_count[t];
      JI = SM / (SA + SB - SM);
      fprintf(outfile, "\t%.10lg", JI);
      break;

    default:
      X = repertoire_matrix[set2_repertoires * s + t];
      fprintf(outfile, "\t%.10lg", X);
      break;
    }
}

uint64_t check_duplicates(struct db * d)
{
  /*
    Test if there are any exact duplicates in the input
    (same repertoire id, same sequence, same V-gene, same J-gene)

    The Zobrist hashing must have already been set up.
  */

  uint64_t sequences = db_getsequencecount(d);
  struct hashtable_s * ht = hash_init(sequences);
  uint64_t dup = 0;

  progress_init("Check duplicates: ", sequences);
  for(uint64_t i=0; i < sequences; i++)
    {
      if (hash_insert(d, ht, nullptr, i))
        dup++;
      progress_update(i);
    }
  progress_done();

  hash_exit(ht);
  ht = nullptr;

  return dup;
}

void overlap(char * set1_filename, char * set2_filename)
{
  /* find overlaps between repertoires */

  db_init();


  /**** Set 1 ****/

  fprintf(logfile, "Immune receptor repertoire set 1\n\n");

  d1 = db_create();
  db_read(d1, set1_filename, opt_existence, "1");

  set1_longestsequence = db_getlongestsequence(d1);
  set1_sequences = db_getsequencecount(d1);
  set1_repertoires = db_get_repertoire_count(d1);
  set1_residues = db_getresiduescount(d1);

  fprintf(logfile, "\n");

  uint64_t set1_sum_size = 0;
  uint64_t set1_sum_count = 0;

  /* determine number of sequences in each of the repertoires (Set 1) */

  set1_repertoire_size = static_cast<uint64_t *>
    (xmalloc(sizeof(uint64_t) * set1_repertoires));
  for (unsigned int s = 0; s < set1_repertoires ; s++)
    set1_repertoire_size[s] = 0;

  set1_repertoire_count = static_cast<uint64_t *>
    (xmalloc(sizeof(uint64_t) * set1_repertoires));
  for (unsigned int s = 0; s < set1_repertoires ; s++)
    set1_repertoire_count[s] = 0;

  set1_repertoire_sq_count = static_cast<double *>
    (xmalloc(sizeof(uint64_t) * set1_repertoires));
  for (unsigned int s = 0; s < set1_repertoires ; s++)
    set1_repertoire_sq_count[s] = 0;

  for (uint64_t i = 0; i < set1_sequences ; i++)
    {
      unsigned int s = db_get_repertoire_id_no(d1, i);
      set1_repertoire_size[s]++;
      uint64_t c = db_get_count(d1, i);
      set1_repertoire_count[s] += c;
      set1_repertoire_sq_count[s] += c * c;
    }

  /* set 1 : sort repertoires alphanumerically for display */

  set1_lookup_repertoire =
    (unsigned int *) xmalloc(sizeof(unsigned int) * set1_repertoires);
  for (unsigned int i = 0; i < set1_repertoires; i++)
    set1_lookup_repertoire[i] = i;

  qsort(set1_lookup_repertoire,
        set1_repertoires,
        sizeof(unsigned int),
        set1_compare_by_repertoire_name);

  /* list of repertoires in set 1 */

  for (unsigned int i = 0; i < set1_repertoires; i++)
    {
      unsigned int s = set1_lookup_repertoire[i];
      set1_sum_size += set1_repertoire_size[s];
      set1_sum_count += set1_repertoire_count[s];
    }

  int w1 = MAX(1, 1 + floor(log10(set1_repertoires)));
  int w2 = MAX(9, 1 + floor(log10(set1_sum_size)));
  int w3 = MAX(5, 1 + floor(log10(set1_sum_count)));

  fprintf(logfile, "Repertoires in set:\n");
  fprintf(logfile, "%*s %*s %*s %s\n",
          w1, "#",
          w2, "Sequences",
          w3, "Count",
          "Repertoire ID");
  for (unsigned int i = 0; i < set1_repertoires; i++)
    {
      unsigned int s = set1_lookup_repertoire[i];
      fprintf(logfile, "%*u %*" PRIu64 " %*" PRIu64 " %s\n",
              w1, i+1,
              w2, set1_repertoire_size[s],
              w3, set1_repertoire_count[s],
              db_get_repertoire_id(d1, s));
    }
  fprintf(logfile, "\n");

  if (opt_existence)
    {
      if (set1_repertoires > 1)
        fatal("Multiple repertoires are not allowed in the first file specified on the command line with the -x or --existence command.");
    }

  /**** Set 2 ****/

  uint64_t set2_sum_size = 0;
  uint64_t set2_sum_count = 0;

  fprintf(logfile, "Immune receptor repertoire set 2\n\n");

  if (set2_filename && strcmp(set1_filename, set2_filename))
    {
      d2 = db_create();
      db_read(d2, set2_filename, false, "2");

      set2_longestsequence = db_getlongestsequence(d2);
      set2_sequences = db_getsequencecount(d2);
      set2_repertoires = db_get_repertoire_count(d2);
      set2_residues = db_getresiduescount(d2);

      fprintf(logfile, "\n");

      /* determine number of sequences in each of the repertoires (Set 2) */

      set2_repertoire_size = static_cast<uint64_t *>
        (xmalloc(sizeof(uint64_t) * set2_repertoires));
      for (unsigned int t = 0; t < set2_repertoires ; t++)
        set2_repertoire_size[t] = 0;

      set2_repertoire_count = static_cast<uint64_t *>
        (xmalloc(sizeof(uint64_t) * set2_repertoires));
      for (unsigned int t = 0; t < set2_repertoires ; t++)
        set2_repertoire_count[t] = 0;

      set2_repertoire_sq_count = static_cast<double *>
        (xmalloc(sizeof(uint64_t) * set2_repertoires));
      for (unsigned int t = 0; t < set2_repertoires ; t++)
        set2_repertoire_sq_count[t] = 0;

      for (uint64_t j = 0; j < set2_sequences ; j++)
        {
          unsigned int t = db_get_repertoire_id_no(d2, j);
          set2_repertoire_size[t]++;
          uint64_t c = db_get_count(d2, j);
          set2_repertoire_count[t] += c;
          set2_repertoire_sq_count[t] += c * c;
        }

      /* set 2 : sort repertoires alphanumerically for display */

      set2_lookup_repertoire =
        (unsigned int *) xmalloc(sizeof(unsigned int) * set2_repertoires);
      for (unsigned int j = 0; j < set2_repertoires; j++)
        set2_lookup_repertoire[j] = j;

      qsort(set2_lookup_repertoire,
            set2_repertoires,
            sizeof(unsigned int),
            set2_compare_by_repertoire_name);

      /* list of repertoires in set 2 */

      for (unsigned int i = 0; i < set2_repertoires; i++)
        {
          unsigned int s = set2_lookup_repertoire[i];
          set2_sum_size += set2_repertoire_size[s];
          set2_sum_count += set2_repertoire_count[s];
        }

      int w1 = MAX(1, 1 + floor(log10(set2_repertoires)));
      int w2 = MAX(9, 1 + floor(log10(set2_sum_size)));
      int w3 = MAX(5, 1 + floor(log10(set2_sum_count)));

      if (set2_repertoires > 0)
        {
          fprintf(logfile, "Repertoires in set:\n");
          fprintf(logfile, "%*s %*s %*s %s\n",
                  w1, "#",
                  w2, "Sequences",
                  w3, "Count",
                  "Repertoire ID");
          for (unsigned int i = 0; i < set2_repertoires; i++)
            {
              unsigned int s = set2_lookup_repertoire[i];
              fprintf(logfile, "%*u %*" PRIu64 " %*" PRIu64 " %s\n",
                      w1, i+1,
                      w2, set2_repertoire_size[s],
                      w3, set2_repertoire_count[s],
                      db_get_repertoire_id(d2, s));
            }
          fprintf(logfile, "\n");
        }
      else
        {
          fatal("Repertoire set missing repertoire_id.");
        }
    }
  else
    {
      /* set2 = set1 */

      d2 = d1;

      fprintf(logfile, "Set 2 is identical to set 1\n");
      fprintf(logfile, "\n");

      set2_longestsequence = db_getlongestsequence(d2);
      set2_sequences = db_getsequencecount(d2);
      set2_repertoires = db_get_repertoire_count(d2);
      set2_residues = db_getresiduescount(d2);

      set2_sum_size = set1_sum_size;
      set2_sum_count = set1_sum_count;

      set2_repertoire_size = set1_repertoire_size;
      set2_repertoire_count = set1_repertoire_count;
      set2_repertoire_sq_count = set1_repertoire_sq_count;
      set2_lookup_repertoire = set1_lookup_repertoire;

      if (set2_repertoires == 0)
        {
          fatal("Repertoire set is missing repertoire_id.");
        }
    }

  unsigned int overall_longest = MAX(set1_longestsequence,
                                     set2_longestsequence);

  fprintf(logfile, "Unique V genes:    %" PRIu64 "\n",
          db_get_v_gene_count());

  fprintf(logfile, "Unique J genes:    %" PRIu64 "\n",
          db_get_j_gene_count());

  /* compute hashes for each sequence in database */

  if (opt_differences <= MAXDIFF_HASH)
    {
      zobrist_init(overall_longest + MAX_INSERTS,
                   db_get_v_gene_count(),
                   db_get_j_gene_count());

      db_hash(d1);

      if (d2 != d1)
        {
          uint64_t dup1 = check_duplicates(d1);
          if (dup1 > 0)
            fprintf(logfile, "Warning: %" PRIu64 " duplicates detected in repertoire set 1\n",
                    dup1);

          db_hash(d2);
        }

      /* store sequences in a hash table */
      /* use an additional bloom filter for increased speed */
      /* hashing into hash table & bloom filter */
      /* check for duplicates in set 2 */

      uint64_t dup2 = 0;
      hashtable = hash_init(set2_sequences);
      bloom_a = bloom_init(hash_get_tablesize(hashtable));
      progress_init("Hashing sequences:", set2_sequences);
      for(uint64_t i=0; i < set2_sequences; i++)
        {
          if (hash_insert(d2, hashtable, bloom_a, i))
            dup2++;
          progress_update(i);
        }
      progress_done();
      if (dup2 > 0)
        fprintf(logfile, "Warning: %" PRIu64 " duplicates detected in repertoire set 2\n", dup2);
    }

  if (! opt_no_matrix)
    {
      if (opt_matrix)
        {
          /* allocate matrix of repertoire set 1 x repertoire set 2 counts */

          repertoire_matrix = static_cast<m_val_t *>
            (xmalloc(sizeof(m_val_t) * set1_repertoires * set2_repertoires));

          for(unsigned int s = 0; s < set1_repertoires; s++)
            for(unsigned int t = 0; t < set2_repertoires; t++)
              repertoire_matrix[set2_repertoires * s + t] = 0;
        }
      else
        {
          /* allocate matrix of set 1 sequences x repertoire set 2 counts */

          repertoire_matrix = static_cast<m_val_t *>
            (xmalloc(sizeof(m_val_t) * set1_sequences * set2_repertoires));

          for(unsigned int s = 0; s < set1_sequences; s++)
            for(unsigned int t = 0; t < set2_repertoires; t++)
              repertoire_matrix[set2_repertoires * s + t] = 0;
        }
    }

  /* compare all sequences */

  pthread_mutex_init(&network_mutex, nullptr);
  pthread_mutex_init(&pairs_mutex, nullptr);
  progress_init("Analysing:        ", set1_sequences);

  if (opt_pairs)
    {
      fprintf(pairsfile,
              "#repertoire_id_1\tsequence_id_1\t"
              "duplicate_count_1\tv_call_1\tj_call_1\t%s_1",
              seq_header);
      for (int k = 0; k < keep_columns_count; k++)
        fprintf(pairsfile, "\t%s_1", keep_columns_names[k]);
      fprintf(pairsfile,
              "\trepertoire_id_2\tsequence_id_2\t"
              "duplicate_count_2\tv_call_2\tj_call_2\t%s_2",
              seq_header);
      for (int k = 0; k < keep_columns_count; k++)
        fprintf(pairsfile, "\t%s_2", keep_columns_names[k]);
      if (opt_distance)
        fprintf(pairsfile, "\tdistance");
      fprintf(pairsfile, "\n");
    }
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
  pthread_mutex_destroy(&pairs_mutex);
  pthread_mutex_destroy(&network_mutex);

  /* dump similarity matrix */

  if (! opt_no_matrix)
    {
      unsigned int x = 0;
      if (opt_alternative)
        {
          if (opt_matrix)
            {
              /* Overlap results, 3-column format */
              progress_init("Writing results:  ",
                            set1_repertoires * set2_repertoires);
              fprintf(outfile, "#repertoire_id_1\trepertoire_id_2\tmatches\n");
              for (unsigned int i = 0; i < set1_repertoires; i++)
                {
                  unsigned int s = set1_lookup_repertoire[i];
                  for (unsigned int j = 0; j < set2_repertoires; j++)
                    {
                      unsigned int t = set2_lookup_repertoire[j];
                      fprintf(outfile,
                              "%s\t%s",
                              db_get_repertoire_id(d1, s),
                              db_get_repertoire_id(d2, t));
                      show_matrix_value(s, t);
                      fprintf(outfile, "\n");
                      progress_update(++x);
                    }
                }
            }
          else
            {
              /* Existence results, 3-column format */
              progress_init("Writing results:  ",
                            set1_sequences * set2_repertoires);
              fprintf(outfile, "#sequence_id_1\trepertoire_id_2\tmatches\n");
              for (unsigned int i = 0; i < set1_sequences; i++)
                {
                  for (unsigned int j = 0; j < set2_repertoires; j++)
                    {
                      unsigned int t = set2_lookup_repertoire[j];
                      fprintf(outfile,
                              "%s\t%s",
                              db_get_sequence_id(d1, i),
                              db_get_repertoire_id(d2, t));
                      show_matrix_value(i, t);
                      fprintf(outfile, "\n");
                      progress_update(++x);
                    }
                }
            }
        }
      else
        {
          if (opt_matrix)
            {
              /* Overlap results, matrix format */
              progress_init("Writing results:  ",
                            set1_repertoires * set2_repertoires);
              fprintf(outfile, "#");
              for (unsigned int j = 0; j < set2_repertoires; j++)
                fprintf(outfile, "\t%s", db_get_repertoire_id(d2, set2_lookup_repertoire[j]));
              fprintf(outfile, "\n");
              for (unsigned int i = 0; i < set1_repertoires; i++)
                {
                  unsigned int s = set1_lookup_repertoire[i];
                  fprintf(outfile, "%s", db_get_repertoire_id(d1, s));
                  for (unsigned int j = 0; j < set2_repertoires; j++)
                    {
                      unsigned int t = set2_lookup_repertoire[j];
                      show_matrix_value(s, t);
                      progress_update(++x);
                    }
                  fprintf(outfile, "\n");
                }
            }
          else
            {
              /* Existence results, matrix format */
              progress_init("Writing results:  ",
                            set1_sequences * set2_repertoires);
              fprintf(outfile, "#");
              for (unsigned int j = 0; j < set2_repertoires; j++)
                fprintf(outfile, "\t%s", db_get_repertoire_id(d2, set2_lookup_repertoire[j]));
              fprintf(outfile, "\n");
              for (unsigned int i = 0; i < set1_sequences; i++)
                {
                  fprintf(outfile, "%s", db_get_sequence_id(d1, i));
                  for (unsigned int j = 0; j < set2_repertoires; j++)
                    {
                      unsigned int t = set2_lookup_repertoire[j];
                      show_matrix_value(i, t);
                      progress_update(++x);
                    }
                  fprintf(outfile, "\n");
                }
            }
        }
    }

  progress_done();
  fprintf(logfile, "\n");

  if (repertoire_matrix)
    xfree(repertoire_matrix);
  repertoire_matrix = nullptr;

  if (opt_differences <= MAXDIFF_HASH)
    {
      bloom_exit(bloom_a);
      hash_exit(hashtable);
      zobrist_exit();
    }

  if (d1 != d2)
    {
      xfree(set2_lookup_repertoire);
      xfree(set2_repertoire_count);
      xfree(set2_repertoire_sq_count);
      xfree(set2_repertoire_size);
      db_free(d2);
    }
  else
    {
      set2_lookup_repertoire = nullptr;
      set2_repertoire_count = nullptr;
      set2_repertoire_sq_count = nullptr;
      set2_repertoire_size = nullptr;
      d2 = nullptr;
    }

  xfree(set1_lookup_repertoire);
  xfree(set1_repertoire_sq_count);
  xfree(set1_repertoire_count);
  xfree(set1_repertoire_size);
  db_free(d1);

  db_exit();
}
