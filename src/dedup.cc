/*
    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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

const uint64_t terminal = -1;
const uint64_t done = -2;

static void report(struct db * d,
                   uint64_t seed,
                   uint64_t * next_seq)
{
  if (next_seq[seed] == done)
    return;

  uint64_t count = opt_ignore_counts ? 1 : db_get_count(d, seed);
  uint64_t link = next_seq[seed];
  next_seq[seed] = done;
  while (link != terminal)
    {
      count += opt_ignore_counts ? 1 : db_get_count(d, link);
      uint64_t temp = next_seq[link];
      next_seq[link] = done;
      link = temp;
    }

  unsigned int seed_rep_id_no = db_get_repertoire_id_no(d, seed);
  const char * seed_rep_id = db_get_repertoire_id(d, seed_rep_id_no);
  const char * seed_v_gene_name = db_get_v_gene_name(d, seed);
  const char * seed_j_gene_name = db_get_j_gene_name(d, seed);

  fprintf(outfile, "%s", seed_rep_id);
  fprintf(outfile, "\t%" PRIu64, count);
  if (! opt_ignore_genes)
    fprintf(outfile, "\t%s\t%s", seed_v_gene_name, seed_j_gene_name);
  fprintf(outfile, "\t");
  db_fprint_sequence(outfile, d, seed);
  fprintf(outfile, "\n");
}


static bool process(struct db * d,
                    hashtable_s * ht,
                    struct bloom_s * b,
                    uint64_t seed,
                    uint64_t * next_seq)
{
  unsigned int seed_rep_id_no = db_get_repertoire_id_no(d, seed);
  uint64_t seed_v_gene = db_get_v_gene(d, seed);
  uint64_t seed_j_gene = db_get_j_gene(d, seed);
  unsigned char * seed_sequence
    = (unsigned char *) db_getsequence(d, seed);
  unsigned int seed_seqlen
    = db_getsequencelen(d, seed);

  uint64_t last = terminal;

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
          unsigned int hit_rep_id_no = db_get_repertoire_id_no(d, hit);

          if (seed_rep_id_no == hit_rep_id_no)
            {
              /* double check that everything matches */
              unsigned int hit_v_gene = db_get_v_gene(d, hit);
              unsigned int hit_j_gene = db_get_j_gene(d, hit);

              if (opt_ignore_genes ||
                  ((seed_v_gene == hit_v_gene) && (seed_j_gene == hit_j_gene)))
                {
                  unsigned char * hit_sequence
                    = (unsigned char *) db_getsequence(d, hit);
                  unsigned int hit_seqlen
                    = db_getsequencelen(d, hit);

                  if ((seed_seqlen == hit_seqlen) &&
                      ! memcmp(seed_sequence, hit_sequence, seed_seqlen))
                    {
                      last = hit;
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

  if (last != terminal)
    {
      next_seq[last] = seed;
      return true;
    }
  else
    return false;
}

void dedup(char * filename)
{
  /* deduplicate a repertoire set */

  db_init();

  struct db * d1 = db_create();

  db_read(d1, filename, false, "1");

  unsigned int longestsequence = db_getlongestsequence(d1);
  uint64_t sequences = db_getsequencecount(d1);

  fprintf(logfile, "Unique V genes:    %" PRIu64 "\n",
          db_get_v_gene_count());

  fprintf(logfile, "Unique J genes:    %" PRIu64 "\n",
          db_get_j_gene_count());


  /* compute hashes for each sequence in database */

  zobrist_init(longestsequence,
               db_get_v_gene_count(),
               db_get_j_gene_count());

  db_hash(d1);

  /* store sequences in a hash table */
  /* use an additional bloom filter for increased speed */
  /* hashing into hash table & bloom filter */

  /* alloc and init array of flags indicating processed sequences */

  uint64_t * next_seq = (uint64_t *) xmalloc(sequences * sizeof(uint64_t));
  for (uint64_t i = 0; i < sequences; i++)
    next_seq[i] = terminal;

  uint64_t dup_seq = 0;

  hashtable_s * hashtable = hash_init(sequences);
  struct bloom_s * bloom = bloom_init(hash_get_tablesize(hashtable));

  fprintf(outfile, "repertoire_id");
  fprintf(outfile, "\tduplicate_count");
  if (! opt_ignore_genes)
    fprintf(outfile, "\tv_call\tj_call");
  fprintf(outfile, "\t%s\n", seq_header);

  progress_init("Deduplicating:    ", sequences);
  for(uint64_t i=0; i < sequences; i++)
    {
      if (process(d1, hashtable, bloom, i, next_seq))
        dup_seq++;
      progress_update(i);
    }
  progress_done();

  fprintf(logfile, "Duplicates merged: %" PRIu64 "\n", dup_seq);

  progress_init("Writing output:   ", sequences);
  for(uint64_t i=0; i < sequences; i++)
    {
      report(d1, i, next_seq);
      progress_update(i);
    }
  progress_done();

  fprintf(logfile, "\n");

  bloom_exit(bloom);
  hash_exit(hashtable);
  hashtable = nullptr;

  xfree(next_seq);
  next_seq = nullptr;

  zobrist_exit();

  db_free(d1);
  db_exit();
}
