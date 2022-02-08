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

const unsigned int no_cluster = UINT_MAX;

static struct iteminfo_s
{
  unsigned int clusterid;
  unsigned int next;
  unsigned int network_start;
  unsigned int network_count;
} * iteminfo = 0;

static struct clusterinfo_s
{
  unsigned int seed;
  unsigned int size;
} * clusterinfo = 0;

static uint64_t clusterinfo_alloc = 0;

static pthread_mutex_t network_mutex;
static unsigned int * network = 0;
static unsigned int network_count = 0;
static unsigned int network_seq = 0;
static uint64_t network_alloc = 0;
static uint64_t seqcount = 0;

static struct db * d;
static struct bloom_s * bloom = 0;
static hashtable_s * hashtable = 0;

static int compare_cluster(const void * a, const void * b)
{
  clusterinfo_s * x = (clusterinfo_s *) a;
  clusterinfo_s * y = (clusterinfo_s *) b;
  if (x->size > y->size)
    return -1;
  else if (x->size < y->size)
    return +1;
  else
    return 0;
}

static inline void hash_insert_cluster(uint64_t seq)
{
  /* find the first empty bucket */
  uint64_t hash = db_gethash(d, seq);
  uint64_t j = hash_getindex(hashtable, hash);
  while (hash_is_occupied(hashtable, j))
    j = hash_getnextindex(hashtable, j);

  hash_set_occupied(hashtable, j);
  hash_set_value(hashtable, j, hash);
  hash_set_data(hashtable, j, seq);
  bloom_set(bloom, hash);
}

static void find_variant_matches(uint64_t seed,
                                 var_s * var,
                                 unsigned int * * hits_data,
                                 unsigned int * hits_count,
                                 uint64_t * hits_alloc)
{
  /* compute hash table index */

  uint64_t j = hash_getindex(hashtable, var->hash);

  /* find matching buckets */

  while (hash_is_occupied(hashtable, j))
    {
      if (hash_compare_value(hashtable, j, var->hash))
        {
          uint64_t hit = hash_get_data(hashtable, j);

          /* double check that everything matches */

          unsigned int seed_v_gene = db_get_v_gene(d, seed);
          unsigned int seed_j_gene = db_get_j_gene(d, seed);

          unsigned int hit_v_gene = db_get_v_gene(d, hit);
          unsigned int hit_j_gene = db_get_j_gene(d, hit);

          if ((seed != hit) &&
              (opt_ignore_genes ||
               ((seed_v_gene == hit_v_gene) &&
                (seed_j_gene == hit_j_gene))))
            {
              unsigned char * seed_sequence
                = (unsigned char *) db_getsequence(d, seed);
              unsigned int seed_seqlen
                = db_getsequencelen(d, seed);
              unsigned char * hit_sequence
                = (unsigned char *) db_getsequence(d, hit);
              unsigned int hit_seqlen
                = db_getsequencelen(d, hit);

              if (check_variant(seed_sequence, seed_seqlen,
                                var,
                                hit_sequence, hit_seqlen))
                {
                  if (*hits_alloc <= *hits_count)
                    {
                      *hits_alloc += 1024;
                      *hits_data = static_cast<unsigned int *>
                        (xrealloc((*hits_data),
                                  (*hits_alloc) * sizeof(unsigned int)));
                    }
                  (*hits_data)[(*hits_count)++] = hit;
                }
            }
        }
      j = hash_getnextindex(hashtable, j);
    }
}

static void process_variants(uint64_t seed,
                             var_s * variant_list,
                             unsigned int * * hits_data,
                             unsigned int * hits_count,
                             uint64_t * hits_alloc)
{
  unsigned int variant_count = 0;
  * hits_count = 0;

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
      if (bloom_get(bloom, var->hash))
        find_variant_matches(seed, var, hits_data, hits_count, hits_alloc);
    }
}

static void process_trad(uint64_t seed,
                         unsigned int * * hits_data,
                         unsigned int * hits_count,
                         uint64_t * hits_alloc)
{
  /* Only to be used with no indels (and d >= 3) */

  for (uint64_t hit = 0; hit < seqcount; hit++)
    if (seed != hit)
      {
        /* check if everything matches */

        unsigned int seed_v_gene = db_get_v_gene(d, seed);
        unsigned int seed_j_gene = db_get_j_gene(d, seed);

        unsigned int hit_v_gene = db_get_v_gene(d, hit);
        unsigned int hit_j_gene = db_get_j_gene(d, hit);

        if (opt_ignore_genes ||
            ((seed_v_gene == hit_v_gene) && (seed_j_gene == hit_j_gene)))
          {
            unsigned int seed_seqlen = db_getsequencelen(d, seed);
            unsigned int hit_seqlen = db_getsequencelen(d, hit);

            if (seed_seqlen == hit_seqlen)
              {
                unsigned char * seed_sequence
                  = (unsigned char *) db_getsequence(d, seed);
                unsigned char * hit_sequence
                  = (unsigned char *) db_getsequence(d, hit);

                if (seq_diff(seed_sequence, hit_sequence, seed_seqlen)
                    <= opt_differences)
                  {
                    if (*hits_alloc <= *hits_count)
                      {
                        *hits_alloc += 1024;
                        *hits_data = static_cast<unsigned int *>
                          (xrealloc((*hits_data),
                                    (*hits_alloc) * sizeof(unsigned int)));
                      }
                    (*hits_data)[(*hits_count)++] = hit;
                  }
              }
          }
      }
}

static void process_seq(uint64_t seed,
                        var_s * variant_list,
                        unsigned int * * hits_data,
                        unsigned int * hits_count,
                        uint64_t * hits_alloc)
{
  if (opt_differences <= MAXDIFF_HASH)
    process_variants(seed, variant_list, hits_data, hits_count, hits_alloc);
  else
    process_trad(seed, hits_data, hits_count, hits_alloc);
}

static void network_thread(int64_t t)
{
  (void) t;

  unsigned int longest = db_getlongestsequence(d);
  uint64_t maxvar = max_variants(longest);

  uint64_t hits_alloc = 1024;
  auto * hits_data = static_cast<unsigned int *>
    (xmalloc(hits_alloc * sizeof(unsigned int)));

  auto * variant_list = static_cast<struct var_s *>
    (xmalloc(maxvar * sizeof(struct var_s)));

  pthread_mutex_lock(&network_mutex);

  while (network_seq < seqcount)
    {
      unsigned int seed = network_seq++;
      progress_update(seed);

      pthread_mutex_unlock(&network_mutex);

      unsigned int hits_count = 0;
      process_seq(seed, variant_list,
                  & hits_data, & hits_count, & hits_alloc);

      pthread_mutex_lock(&network_mutex);

      iteminfo[seed].network_start = network_count;
      iteminfo[seed].network_count = hits_count;

      if (network_count + hits_count > network_alloc)
        {
          while (network_count + hits_count > network_alloc)
            network_alloc += 1024 * 1024;

          network = static_cast<unsigned int*>
            (xrealloc(network, network_alloc * sizeof(unsigned int)));
        }

      for(unsigned int k = 0; k < hits_count; k++)
        network[network_count++] = hits_data[k];
    }

  pthread_mutex_unlock(&network_mutex);

  xfree(variant_list);
  xfree(hits_data);
}

static unsigned int clustersize = 0;
static unsigned int current_cluster_tail = 0;

static void process_seed(unsigned int seed)
{
  clustersize++;

  unsigned int s = iteminfo[seed].network_start;
  unsigned int c = iteminfo[seed].network_count;

  unsigned int clusterid = iteminfo[seed].clusterid;

  for(unsigned int i = 0; i < c; i++)
    {
      unsigned int hit = network[s + i];
      if (iteminfo[hit].clusterid == no_cluster)
        {
          /* add hit to cluster, update linked chain */
          iteminfo[hit].clusterid = clusterid;
          iteminfo[current_cluster_tail].next = hit;
          current_cluster_tail = hit;
        }
    }
}

void cluster(char * filename)
{
  fprintf(logfile, "Immune receptor repertoire clustering\n\n");

  db_init();

  d = db_create();
  db_read(d, filename, false, false);

  unsigned int longest = db_getlongestsequence(d);
  seqcount = db_getsequencecount(d);

  fprintf(logfile, "\n");
  fprintf(logfile, "Unique V genes:    %" PRIu64 "\n",
          db_get_v_gene_count());
  fprintf(logfile, "Unique J genes:    %" PRIu64 "\n",
          db_get_j_gene_count());
  fprintf(logfile, "\n");

  if (opt_differences <= MAXDIFF_HASH)
    {
      zobrist_init(longest + MAX_INSERTS,
                   db_get_v_gene_count(),
                   db_get_j_gene_count());

      db_hash(d);

      hashtable = hash_init(seqcount);
      bloom = bloom_init(hash_get_tablesize(hashtable) * 2);
    }

  iteminfo = static_cast<struct iteminfo_s *>
    (xmalloc(seqcount * sizeof(struct iteminfo_s)));

  progress_init("Hashing sequences:", seqcount);
  for(uint64_t i=0; i < seqcount; i++)
    {
      iteminfo[i].clusterid = no_cluster;
      iteminfo[i].next = no_cluster;
      if (opt_differences <= MAXDIFF_HASH)
        hash_insert_cluster(i);
      progress_update(i);
    }
  progress_done();

  network = static_cast<unsigned int*>
    (xmalloc(network_alloc * sizeof(unsigned int)));
  network_count = 0;
  network_seq = 0;

  pthread_mutex_init(&network_mutex, nullptr);
  progress_init("Building network: ", seqcount);

  if (opt_threads == 1)
    {
      network_thread(0);
    }
  else
    {
      ThreadRunner * sim_tr = new ThreadRunner(static_cast<int>(opt_threads),
                                               network_thread);
      sim_tr->run();
      delete sim_tr;
    }

  progress_done();
  pthread_mutex_destroy(&network_mutex);


  unsigned int clustercount = 0;

  progress_init("Clustering:       ", seqcount);

  /* for each non-clustered item, look for subseeds ... */
  uint64_t x = 0;
  for(unsigned int seed = 0; seed < seqcount; seed++)
    {
      struct iteminfo_s * ap = iteminfo + seed;

      if (ap->clusterid == no_cluster)
        {
          /* start a new cluster with a new initial seed */

          ap->clusterid = clustercount;
          ap->next = no_cluster;
          current_cluster_tail = seed;
          clustersize = 0;

          /* find initial matches */
          process_seed(seed);
          progress_update(++x);

          unsigned int subseed = ap->next;

          /* process all subseeds */
          while(subseed != no_cluster)
            {
              process_seed(subseed);
              progress_update(++x);
              subseed = iteminfo[subseed].next;
            }

          if (clustercount >= clusterinfo_alloc)
            {
              /* allocate memory for more clusters... */
              clusterinfo_alloc += 1024;
              clusterinfo = static_cast<struct clusterinfo_s *>
                (xrealloc(clusterinfo,
                          clusterinfo_alloc * sizeof(clusterinfo_s)));
            }

          struct clusterinfo_s * sp = clusterinfo + clustercount;
          sp->seed = seed;
          sp->size = clustersize;
          clustercount++;
        }
    }

  progress_done();

  progress_init("Sorting clusters: ", clustercount);
  qsort(clusterinfo, clustercount, sizeof(clusterinfo_s), compare_cluster);
  progress_done();

  /* dump clusters */

  uint64_t j = 0;
  progress_init("Writing clusters: ", seqcount);
  fprintf(outfile,
          "#cluster_no\tcluster_size\trepertoire_id\tsequence_id\t"
          "duplicate_count\tv_call\tj_call\t%s\n",
          opt_nucleotides ? "junction" : "junction_aa");
  for(unsigned int i = 0; i < clustercount; i++)
    {
      unsigned int seed = clusterinfo[i].seed;
      unsigned int size = clusterinfo[i].size;
      for(unsigned int a = seed; a != no_cluster; a = iteminfo[a].next)
        {
          fprintf(outfile,
                  "%u\t%u\t",
                  i + 1,
                  size);
          fprintf(outfile,
                  "%s\t%s\t%" PRIu64 "\t%s\t%s\t",
                  db_get_repertoire_id(d, db_get_repertoire_id_no(d, a)),
                  db_get_sequence_id(d, a),
                  db_get_count(d, a),
                  db_get_v_gene_name(d, a),
                  db_get_j_gene_name(d, a));
          db_fprint_sequence(outfile, d, a);
          fprintf(outfile, "\n");
          j++;
        }
      progress_update(j);
    }
  progress_done();

  fprintf(logfile, "\n");
  fprintf(logfile, "Clusters:          %u\n", clustercount);

  xfree(network);
  if (clusterinfo)
    xfree(clusterinfo);
  if (iteminfo)
    xfree(iteminfo);

  if (opt_differences <= MAXDIFF_HASH)
    {
      bloom_exit(bloom);
      hash_exit(hashtable);
      zobrist_exit();
    }

  db_free(d);
  db_exit();
}
