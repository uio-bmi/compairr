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

//#define DUMP 1
#define SMART 1
//#define USEBLOOM 1
//#define USEHASH 1
//#define HASHSTATS 1

#ifdef HASHSTATS
static uint64_t tries = 0;
#endif

#ifdef DUMP
static unsigned char * dump_seed_sequence;
static unsigned int dump_seed_seqlen;
static uint64_t dump_v_gene;
static uint64_t dump_d_gene;
#endif

void ps(unsigned int seqlen, unsigned char * sequence)
{
  char ab[] = "ACDEFGHIKLMNPQRSTVWY";
  for(unsigned int i=0; i < seqlen; i++)
    printf("%c", ab[sequence[i]]);
}

inline void seq_copy(unsigned char * a,
                     unsigned int a_start,
                     unsigned char * b,
                     unsigned int b_start,
                     unsigned int length)
{
  /* copy part of the sequence b to a */
  memcpy(a + a_start, b + b_start, length);
}

inline bool seq_identical(unsigned char * a,
                          unsigned int a_start,
                          unsigned char * b,
                          unsigned int b_start,
                          unsigned int length)
{
  /* compare parts of two sequences a and b */
  /* return false if different, true if identical */

  return ! memcmp(a + a_start, b + b_start, length);
}

uint64_t max_variants(uint64_t longest)
{
  /*
    calculate upper limit on number of possible variants

    d=0
    identical:     1

    d=1
    deletions:     longest
    substitutions: (ab_size - 1) * longest
    insertions:    (ab_size - 1) * (longest + 1) + 1
  */

  uint64_t maxvar = 1; // identical non-variant

  for (int i = 0; i < opt_differences; i++)
    maxvar *= (2 * (longest + i) + 1) * alphabet_size - (longest + i);

  return maxvar;
}

void generate_variant_sequence(unsigned char * seed_sequence,
                               unsigned int seed_seqlen,
                               struct var_s * var,
                               unsigned char * seq,
                               unsigned int * seqlen)
{
  /* generate the actual sequence of a variant */

  switch (var->kind)
    {
    case identical:
      seq_copy(seq, 0, seed_sequence, 0, seed_seqlen);
      * seqlen = seed_seqlen;
      break;

    case substitution:
      seq_copy(seq, 0, seed_sequence, 0, seed_seqlen);
      seq[var->pos1] = var->residue1;
      * seqlen = seed_seqlen;
      break;

    case deletion:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos1);
      seq_copy(seq, var->pos1,
               seed_sequence, var->pos1 + 1,
               seed_seqlen - var->pos1 - 1);
      * seqlen = seed_seqlen - 1;
      break;

    case insertion:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos1);
      seq[var->pos1] = var->residue1;
      seq_copy(seq, var->pos1 + 1,
               seed_sequence, var->pos1,
               seed_seqlen - var->pos1);
      * seqlen = seed_seqlen + 1;
      break;

    case sub_sub:
      seq_copy(seq, 0, seed_sequence, 0, seed_seqlen);
      seq[var->pos1] = var->residue1;
      seq[var->pos2] = var->residue2;
      * seqlen = seed_seqlen;
      break;

    case del_del:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos1);
      seq_copy(seq, var->pos1,
               seed_sequence, var->pos1 + 1,
               var->pos2 - var->pos1 - 1);
      seq_copy(seq, var->pos2 - 1,
               seed_sequence, var->pos2 + 1,
               seed_seqlen - var->pos2 - 1);
      * seqlen = seed_seqlen - 2;
      break;

    case ins_ins:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos1);
      seq[var->pos1] = var->residue1;
      seq_copy(seq, var->pos1 + 1,
               seed_sequence, var->pos1,
               var->pos2 - var->pos1 - 1);
      seq[var->pos2] = var->residue2;
      seq_copy(seq, var->pos2 + 1,
               seed_sequence, var->pos2 - 1,
               seed_seqlen - var->pos2 + 1);
      * seqlen = seed_seqlen + 2;
      break;

    case ins_sub:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos1);
      seq_copy(seq, var->pos1 + 1,
               seed_sequence, var->pos1,
               seed_seqlen - var->pos1);
      seq[var->pos1] = var->residue1;
      seq[var->pos2] = var->residue2;
      * seqlen = seed_seqlen + 1;
      break;

    case del_sub:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos1);
      seq_copy(seq, var->pos1,
               seed_sequence, var->pos1 + 1,
               seed_seqlen - var->pos1 - 1);
      seq[var->pos2] = var->residue2;
      * seqlen = seed_seqlen - 1;
      break;

    case del_ins:
      if (var->pos1 <= var->pos2)
        {
          // deletion left of insertion
          seq_copy(seq, 0,
                   seed_sequence, 0,
                   var->pos1);
          seq_copy(seq, var->pos1,
                   seed_sequence, var->pos1 + 1,
                   var->pos2 - var->pos1);
          seq_copy(seq, var->pos2 + 1,
                   seed_sequence, var->pos2 + 1,
                   seed_seqlen - var->pos2 - 1);
          seq[var->pos2] = var->residue2;
        }
      else
        {
          // insertion left of deletion
          seq_copy(seq, 0,
                   seed_sequence, 0,
                   var->pos2);
          seq_copy(seq, var->pos2 + 1,
                   seed_sequence, var->pos2,
                   var->pos1 - var->pos2);
          seq_copy(seq, var->pos1 + 1,
                   seed_sequence, var->pos1 + 1,
                   seed_seqlen - var->pos1 - 1);
          seq[var->pos2] = var->residue2;
        }
      * seqlen = seed_seqlen;
      break;

    case sub_del:
    case sub_ins:
    case ins_del:
    default:
      fatal("Internal error");
      break;

    }
}


bool check_variant(unsigned char * seed_sequence,
                   unsigned int seed_seqlen,
                   var_s * var,
                   unsigned char * amp_sequence,
                   unsigned int amp_seqlen)
{
  return true; // just ignore the possibility

  /* make sure seed with given variant is really identical to amp */
  /* we know the hashes are identical */

  if (var->kind > 4)
    fprintf(stderr, "Var kind %d\n", var->kind);

  bool equal = false;

  switch (var->kind)
    {
    case identical:
      equal = ((seed_seqlen == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              seed_seqlen)));
      break;

    case substitution:
      equal = ((seed_seqlen == amp_seqlen) &&
               (amp_sequence[var->pos1] == var->residue1) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos1)) &&
               (seq_identical(seed_sequence, var->pos1 + 1,
                              amp_sequence,  var->pos1 + 1,
                              seed_seqlen - var->pos1 - 1)));
      break;

    case deletion:
      equal = (((seed_seqlen - 1) == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos1)) &&
               (seq_identical(seed_sequence, var->pos1 + 1,
                              amp_sequence,  var->pos1,
                              seed_seqlen - var->pos1 - 1)));
      break;

    case insertion:
      equal = (((seed_seqlen + 1) == amp_seqlen) &&
               (amp_sequence[var->pos1] == var->residue1) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos1)) &&
               (seq_identical(seed_sequence, var->pos1,
                              amp_sequence,  var->pos1 + 1,
                              seed_seqlen - var->pos1)));
      break;

    case sub_sub:
      equal = ((seed_seqlen == amp_seqlen) &&
               (amp_sequence[var->pos1] == var->residue1) &&
               (amp_sequence[var->pos2] == var->residue2) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos1)) &&
               (seq_identical(seed_sequence, var->pos1 + 1,
                              amp_sequence, var->pos1 + 1,
                              var->pos2 - var->pos1 - 1)) &&
               (seq_identical(seed_sequence, var->pos2 + 1,
                              amp_sequence,  var->pos2 + 1,
                              seed_seqlen - var->pos2 - 1)));
      break;

    case del_del:
    case ins_ins:
    case sub_del:
    case sub_ins:
    case del_sub:
    case ins_sub:
    case ins_del:
    case del_ins:
    default:
      fatal("Internal error");
      break;

    }

  return equal;
}

inline void add_variant(uint64_t hash,
                        var_s * variant_list,
                        unsigned int * variant_count,
                        enum mutation_kind_enum kind,
                        unsigned int pos1,
                        unsigned char residue1,
                        unsigned int pos2,
                        unsigned char residue2,
                        struct bloom_s * bloom,
                        struct hashtable_s * ht)
{
  (void) bloom;

#ifdef HASHSTATS
  tries++;
#endif

#ifdef USEBLOOM
  if (bloom_get(bloom, hash))
    printf("BLOOM HIT!\n");
#endif

  // check hash table

#ifdef USEHASH
  uint64_t j = hash_getindex(ht, hash);
  while (hash_is_occupied(ht, j))
    {
      if (hash_compare_value(ht, j, hash))
        return;
      j = hash_getnextindex(ht, j);
    }

  hash_set_occupied(ht, j);
  hash_set_value(ht, j, hash);
#else
  (void) ht;
#endif

  var_s * v = variant_list + (*variant_count)++;
  v->hash = hash;
  v->kind = kind;
  v->pos1 = pos1;
  v->residue1 = residue1;
  v->pos2 = pos2;
  v->residue2 = residue2;

#ifdef DUMP

#ifdef SMART
  char kindletters[] = "= s d i sssdsidsdddiisidii";
  unsigned char * mod = (unsigned char *) xmalloc(32);
  unsigned int modlen = 0;
  generate_variant_sequence(dump_seed_sequence,
                            dump_seed_seqlen,
                            v,
                            mod,
                            & modlen);

  uint64_t newhash = zobrist_hash(mod, modlen, dump_v_gene, dump_d_gene);

  printf("calc: %016llx correct: %016llx", hash, newhash);
  printf(" %.2s", kindletters + 2 * kind);
  printf(" pos: %2d,%2d res: ", pos1, pos2);
  ps(1, &residue1);
  ps(1, &residue2);
  printf(" ");
  ps(modlen, mod);
  if (hash != newhash)
    printf(" !!!");
  printf("\n");

  xfree(mod);

#else

  printf("calc: %016llx\n", hash);

#endif

#endif
}

void generate_variants_0(uint64_t hash,
                         unsigned char * sequence,
                         unsigned int seqlen,
                         uint64_t v_gene,
                         uint64_t d_gene,
                         var_s * variant_list,
                         unsigned int * variant_count,
                         struct bloom_s * bloom,
                         struct hashtable_s * ht)
{
  (void) sequence;
  (void) seqlen;
  (void) v_gene;
  (void) d_gene;

  /* identical non-variant */
  add_variant(hash,
              variant_list, variant_count,
              identical, 0, 0, 0, 0,
              bloom, ht);
}

void generate_variants_1(uint64_t hash,
                         unsigned char * sequence,
                         unsigned int seqlen,
                         uint64_t v_gene,
                         uint64_t d_gene,
                         var_s * variant_list,
                         unsigned int * variant_count,
                         struct bloom_s * bloom,
                         struct hashtable_s * ht)
{

#if 1

  /* substitutions */

  for(unsigned int i = 0; i < seqlen; i++)
    {
      unsigned char residue1 = sequence[i];
      uint64_t hash1 = hash ^ zobrist_value(i, residue1);
      for (unsigned char v = 0; v < alphabet_size; v++)
        if (v != residue1)
          {
            uint64_t hash2 = hash1 ^ zobrist_value(i, v);

#ifdef DUMP
            printf("subst %d %d\n", i, v);
#endif

            add_variant(hash2,
                        variant_list, variant_count,
                        substitution, i, v, 0, 0,
                        bloom, ht);
          }
    }

#endif

  /* indels */

  if (opt_indels)
    {

#if 1

      /* deletions */

#ifdef DUMP
      printf("1del\n");
#endif

      hash = zobrist_hash_delete_first
        (reinterpret_cast<unsigned char *> (sequence),
         seqlen,
         v_gene,
         d_gene);
      add_variant(hash,
                  variant_list, variant_count,
                  deletion, 0, 0, 0, 0,
                  bloom, ht);
      unsigned char deleted = sequence[0];
      for(unsigned int i = 1; i < seqlen; i++)
        {
          unsigned char v = sequence[i];
          if (v != deleted)
            {
              hash ^= zobrist_value(i - 1, deleted)
                ^ zobrist_value(i - 1, v);
              add_variant(hash,
                          variant_list, variant_count,
                          deletion, i, 0, 0, 0,
                          bloom, ht);
              deleted = v;
            }
        }

#endif

#if 1

      /* insertions */

#ifdef DUMP
      printf("1ins\n");
#endif

      hash = zobrist_hash_insert_first
        (reinterpret_cast<unsigned char *>(sequence),
         seqlen,
         v_gene,
         d_gene);
      for (unsigned char v = 0; v < alphabet_size; v++)
        {
          uint64_t hash1 = hash ^ zobrist_value(0, v);
          add_variant(hash1,
                      variant_list, variant_count,
                      insertion, 0, v, 0, 0,
                      bloom, ht);
        }
      for (unsigned int i = 0; i < seqlen; i++)
        {
          unsigned char inserted = sequence[i];
          hash ^= zobrist_value(i, inserted) ^ zobrist_value(i+1, inserted);
          for (unsigned char v = 0; v < alphabet_size; v++)
            if (v != inserted)
              {
                uint64_t hash1 = hash ^ zobrist_value(i + 1, v);
                add_variant(hash1,
                            variant_list, variant_count,
                            insertion, i + 1, v, 0, 0,
                            bloom, ht);
              }
        }

#endif

    }
}

void generate_variants_2_all(uint64_t hash,
                             unsigned char * sequence,
                             unsigned int seqlen,
                             uint64_t v_gene,
                             uint64_t d_gene,
                             var_s * variant_list,
                             unsigned int * variant_count,
                             struct bloom_s * bloom,
                             struct hashtable_s * ht)
{
  /* make a copy of the sequence */

  unsigned char * mut = (unsigned char *) xmalloc(seqlen + 1);
  unsigned int mutlen = seqlen;
  memcpy(mut, sequence, mutlen);

  /* substitutions */

  for(unsigned int i = 0; i < seqlen; i++)
    {
      unsigned char residue1 = sequence[i];
      uint64_t hash1 = hash ^ zobrist_value(i, residue1);
      for (unsigned char v = 0; v < alphabet_size; v++)
        if (v != residue1)
          {
            uint64_t hash2 = hash1 ^ zobrist_value(i, v);
            add_variant(hash2,
                        variant_list, variant_count,
                        substitution, i, v, 0, 0,
                        bloom, ht);
            mut[i] = v;

#ifdef DUMP
            printf("sub %2d %2d ", i, v);
            ps(mutlen, mut);
            printf("\n");
#endif

            generate_variants_1(hash2,
                                mut, mutlen,
                                v_gene, d_gene,
                                variant_list, variant_count,
                                bloom, ht);
          }
      mut[i] = residue1;
    }

  /* indels */

  if (opt_indels)
    {
      /* deletions */

      mutlen = seqlen - 1;
      memcpy(mut, sequence + 1, mutlen);

      hash = zobrist_hash_delete_first
        (reinterpret_cast<unsigned char *> (sequence),
         seqlen,
         v_gene,
         d_gene);

      add_variant(hash,
                  variant_list, variant_count,
                  deletion, 0, 0, 0, 0,
                  bloom, ht);

#ifdef DUMP
      printf("del %2d ", 0);
      ps(mutlen, mut);
      printf("\n");
#endif

      generate_variants_1(hash,
                          mut, mutlen,
                          v_gene, d_gene,
                          variant_list, variant_count,
                          bloom, ht);

      unsigned char deleted = sequence[0];
      for(unsigned int i = 1; i < seqlen; i++)
        {
          unsigned char v = sequence[i];
          if (v != deleted)
            {
              hash ^= zobrist_value(i - 1, deleted)
                ^ zobrist_value(i - 1, v);
              add_variant(hash,
                          variant_list, variant_count,
                          deletion, i, 0, 0, 0,
                          bloom, ht);
              mut[i - 1] = sequence[i-1];

#ifdef DUMP
              printf("del %2d ", i);
              ps(mutlen, mut);
              printf("\n");
#endif

              generate_variants_1(hash,
                                  mut, mutlen,
                                  v_gene, d_gene,
                                  variant_list, variant_count,
                                  bloom, ht);
              mut[i - 1] = deleted;
              deleted = v;
            }
        }

#if 1
      /* insertions */

      mutlen = seqlen + 1;
      memcpy(mut + 1, sequence, mutlen - 1);

      hash = zobrist_hash_insert_first
        (reinterpret_cast<unsigned char *>(sequence),
         seqlen,
         v_gene,
         d_gene);
      for (unsigned char v = 0; v < alphabet_size; v++)
        {
          uint64_t hash1 = hash ^ zobrist_value(0, v);
          add_variant(hash1,
                      variant_list, variant_count,
                      insertion, 0, v, 0, 0,
                      bloom, ht);
          mut[0] = v;

#ifdef DUMP
          printf("ins %2d %2d ", 0, v);
          ps(mutlen, mut);
          printf("\n");
#endif

          generate_variants_1(hash1,
                              mut, mutlen,
                              v_gene, d_gene,
                              variant_list, variant_count,
                              bloom, ht);
        }

      mut[0] = sequence[0];

      for (unsigned int i = 0; i < seqlen; i++)
        {
          unsigned char inserted = sequence[i];
          hash ^= zobrist_value(i, inserted) ^ zobrist_value(i+1, inserted);
          for (unsigned char v = 0; v < alphabet_size; v++)
            if (v != inserted)
              {
                uint64_t hash1 = hash ^ zobrist_value(i + 1, v);
                add_variant(hash1,
                            variant_list, variant_count,
                            insertion, i + 1, v, 0, 0,
                            bloom, ht);
                mut[i + 1] = v;

#ifdef DUMP
                printf("ins %2d %2d ", i+1, v);
                ps(mutlen, mut);
                printf("\n");
#endif

                generate_variants_1(hash1,
                                    mut, mutlen,
                                    v_gene, d_gene,
                                    variant_list, variant_count,
                                    bloom, ht);
              }
          mut[i + 1] = sequence[i + 1];
        }
#endif
    }

  xfree(mut);
}

void generate_variants_2_smart(uint64_t hash,
                               unsigned char * sequence,
                               unsigned int seqlen,
                               uint64_t v_gene,
                               uint64_t d_gene,
                               var_s * variant_list,
                               unsigned int * variant_count,
                               struct bloom_s * bloom,
                               struct hashtable_s * ht)
{
  (void) v_gene;
  (void) d_gene;

  /* generate all double substitutions */

#if 1

  for (unsigned int i = 0; i < seqlen; i++)
    {
      unsigned char res1 = sequence[i];
      uint64_t hash1 = hash ^ zobrist_value(i, res1);

      for (unsigned char v = 0; v < alphabet_size; v++)
        {
          if (v != res1)
            {
              uint64_t hash2 = hash1 ^ zobrist_value(i, v);

              for (unsigned int j = i + 1; j < seqlen; j++)
                {
                  unsigned char res2 = sequence[j];
                  uint64_t hash3 = hash2 ^ zobrist_value(j, res2);

                  for (unsigned char w = 0; w < alphabet_size; w++)
                    {
                      if (w != res2)
                        {
                          uint64_t hash4 = hash3 ^ zobrist_value(j, w);
                          add_variant(hash4,
                                      variant_list, variant_count,
                                      sub_sub, i, v, j, w,
                                      bloom, ht);
                        }
                    }
                }
            }
        }
    }

#endif

  if (opt_indels)
    {

#if 1

      /* generate all double deletions */
      /* avoid special cases that result in identical sequences */
      /* e.g. with repeats like AAC => --C = C, -A- = A, A-- = A */
      /* or like ATAT => --AT = AT, -T-T = TT, -TA- = TA, A--T = AT,
                         A-A- = AA, AT-- = AT */

      if (seqlen >= 2)
        {
          uint64_t hash1 = zobrist_hash_delete_first_two
            (reinterpret_cast<unsigned char *> (sequence),
             seqlen,
             v_gene,
             d_gene);

          for (unsigned int i = 0; i < seqlen - 1; i++)
            {
              if (i > 0)
                {
                  hash1 ^= zobrist_value(i - 1, sequence[i - 1])
                    ^ zobrist_value(i - 1, sequence[i + 1]);
                }

              if ((i == 0) ||
                  ((sequence[i - 1] != sequence[i + 1]) &&
                   (sequence[i - 1] != sequence[i    ])))
                {
                  uint64_t hash2 = hash1;

                  for (unsigned int j = i + 1; j < seqlen; j++)
                    {
                      if (j > i + 1)
                        {
                          hash2 ^= zobrist_value(j - 2, sequence[j - 1])
                            ^ zobrist_value(j - 2, sequence[j]);
                        }

                      if ((j <= i + 1) || (sequence[j - 1] != sequence[j]))
                        {
                          add_variant(hash2,
                                      variant_list, variant_count,
                                      del_del, i, 0, j, 0,
                                      bloom, ht);
                        }
                    }
                }
            }
        }
#endif

      uint64_t hash1;

#if 1

      /* generate all double insertions */

      hash1 = zobrist_hash_insert_first_two
        (reinterpret_cast<unsigned char *> (sequence),
         seqlen,
         v_gene,
         d_gene);

      /* insert any residue in first position */
      /* except from the same as in the previous position */

      for (unsigned int i = 0; i < seqlen + 1; i++)
        {
          if (i > 0)
            {
              hash1 ^= zobrist_value(i - 1, sequence[i - 1]) ^
                       zobrist_value(i + 1, sequence[i - 1]);
            }
          for (unsigned char v1 = 0; v1 < alphabet_size; v1++)
            {
              if ((i == 0) || (v1 != sequence[i - 1]))
                {
                  uint64_t hash2 = hash1 ^ zobrist_value(i, v1);

                  /* insert any residue in second position */
                  /* except from the same as in the previous position */

                  for (unsigned int j = i + 1; j < seqlen + 2; j++)
                    {
                      if (j > i + 1)
                        {
                          hash2 ^= zobrist_value(j - 1, sequence[j - 2]) ^
                            zobrist_value(j,     sequence[j - 2]);
                        }

                      for (unsigned char v2 = 0; v2 < alphabet_size; v2++)
                        {
                          if ((j == i + 1) || (v2 != sequence[j - 2]))
                            {
                              uint64_t hash3 = hash2 ^ zobrist_value(j, v2);
                              add_variant(hash3,
                                          variant_list, variant_count,
                                          ins_ins, i, v1, j, v2,
                                          bloom, ht);
                            }
                        }
                    }
                }
            }
        }

#endif


#if 1

      /* generate all insertions/substitutions */

      hash1 = zobrist_hash_insert_first
        (reinterpret_cast<unsigned char *> (sequence),
         seqlen,
         v_gene,
         d_gene);

      for (unsigned int i = 0; i < seqlen + 1; i++)
        {
          if (i > 0)
            {
              hash1 ^= zobrist_value(i - 1, sequence[i - 1]) ^
                       zobrist_value(i,     sequence[i - 1]);
            }
          for (unsigned char v1 = 0; v1 < alphabet_size; v1++)
            {
              if (((i == 0) )// && (sequence[0] != sequence[1]))
                  || ((i > 0) && (v1 != sequence[i - 1])))
                {
                  /* insert any residue in position i */
                  /* except the same as in the previous position */

                  uint64_t hash2 = hash1 ^ zobrist_value(i, v1);

                  for (unsigned int j = 0; j < seqlen + 1; j++)
                    {
                      if ((j < i) || (j > i + 1))
                        {
                          /*
                            Substitute any residue in position j.
                            Avoid same position as insertion,
                            as it would be same as single insertion.
                            Also avoid position right after insertion,
                            as it would give the same result.
                          */

                          /* p is the position of the substitution
                             in the original sequence */

                          unsigned int p = j;
                          if (i < j)
                            p--;

                          if ((i + 1 != j) || (sequence[i] != sequence[j]))
                            {
                              for (unsigned char v2 = 0;
                                   v2 < alphabet_size;
                                   v2++)
                                {
                                  /* avoid same residue as original */
                                  if (v2 != sequence[p])
                                    {
                                      uint64_t hash3 = hash2
                                        ^ zobrist_value(j, v2)
                                        ^ zobrist_value(j, sequence[p]);
                                      add_variant(hash3,
                                                  variant_list, variant_count,
                                                  ins_sub, i, v1, j, v2,
                                                  bloom, ht);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

#endif

#if 1

      /* generate all deletions/substitutions */
      /* still not 100% correct, e.g. with sequence CCAC */

      hash1 = zobrist_hash_delete_first
        (reinterpret_cast<unsigned char *> (sequence),
         seqlen,
         v_gene,
         d_gene);

      for (unsigned int i = 0; i < seqlen; i++)
        {
          /* i = zero-based position of deletion in orig sequence */
          unsigned char deleted = sequence[i];
          unsigned int rep = 0;

          if ((i == 0) || (sequence[i-1] != sequence[i]))
            {
              if (i > 0)
                {
                  hash1 ^= zobrist_value(i - 1, sequence[i - 1])
                         ^ zobrist_value(i - 1, sequence[i    ]);
                }

              for (unsigned int j = 0; j < seqlen; j++)
                {
                  if (j == i + 3)
                    rep = 1;
                  if (rep && (sequence[j-1] != sequence[j-2]))
                    rep = 0;

                  /* j = zero-based position of subst in orig sequence */

                  if ((i > j + 1) || (i < j))
                    {
                      /*
                        Substitute any residue in position j.
                        Avoid same position as deletion.
                        Also avoid position left of deletion.
                      */

                      /* p = zero-based position of subst in new sequence */

                      unsigned int p = j;
                      if (i < j)
                        p--;

                      for (unsigned char v2 = 0; v2 < alphabet_size; v2++)
                        {
                          /*
                            When substituting in pos j, avoid same residue as
                            - the original residue in position j
                            - just deleted in position j-1
                            - between deletion and substitution if j=i+2
                            - repeated between i and j when i < j
                          */

                          if ((v2 != sequence[j]) &&
                              ((i + 1 != j) || (v2 != deleted)) &&
                              ((i + 2 != j) || (v2 != sequence[i+1])) &&
                              ((! rep) || (v2 != sequence[j-1])))
                            {
                              uint64_t hash2 = hash1
                                ^ zobrist_value(p, v2)
                                ^ zobrist_value(p, sequence[j]);
                              add_variant(hash2,
                                          variant_list, variant_count,
                                          del_sub, i, 0, p, v2,
                                          bloom, ht);
                            }
                        }
                    }
                }
            }
        }

#endif

#if 1

      /* generate all deletion/insertions */
      /* avoid same as single or double substitutions */

      hash1 = hash ^ zobrist_value(0, sequence[0]);

      for (unsigned int i = 0; i < seqlen; i++)
        {
          /* i is position of deletion in original sequence */

          if (i > 0)
            {
              hash1 ^= zobrist_value(i, sequence[i])
                     ^ zobrist_value(i, sequence[i-1]);
            }

          uint64_t hash2 = hash1;

          for (unsigned int j = 0; j < seqlen; j++)
            {
              /* j is position of insertion in new sequence */

              /* p is position left of insertion in old sequence */

              unsigned int p = j;

              if (j > 0)
                {
                  if (j <= i)
                    p--;
                  hash2 ^= zobrist_value(j,   sequence[p])
                         ^ zobrist_value(j-1, sequence[p]);
                }

              if (j != i)
                {
                  for (unsigned char v = 0; v < alphabet_size; v++)
                    {
                      if ((v != sequence[j]) &&
                          ((j + 1 != i) || (v != sequence[i])))
                        {
                          uint64_t hash3 = hash2 ^ zobrist_value(j, v);
                          add_variant(hash3,
                                      variant_list, variant_count,
                                      del_ins, i, 0, j, v,
                                      bloom, ht);
                        }
                    }
                }
            }
        }

#endif

    }
}

void generate_variants(uint64_t hash,
                       unsigned char * sequence,
                       unsigned int seqlen,
                       uint64_t v_gene,
                       uint64_t d_gene,
                       var_s * variant_list,
                       unsigned int * variant_count)
{
#ifdef DUMP
  dump_seed_sequence = sequence;
  dump_seed_seqlen = seqlen;
  dump_v_gene = v_gene;
  dump_d_gene = d_gene;

#if 0
  printf("seq (%i): ", seqlen);
  ps(seqlen, sequence);
  printf("\n");
#endif
#endif


#ifdef USEBLOOM
  uint64_t bloomsize = 1;
  uint64_t needed = max_variants(seqlen) * 2;
  while (needed > bloomsize)
    bloomsize *= 2;

  struct bloom_s * bloom = bloom_init(bloomsize);
#else
  struct bloom_s * bloom = nullptr;
#endif

#ifdef USEHASH
  struct hashtable_s * ht = hash_init(max_variants(seqlen));
#else
  struct hashtable_s * ht = nullptr;
#endif

#if 1
  generate_variants_0(hash,
                      sequence, seqlen,
                      v_gene, d_gene,
                      variant_list, variant_count,
                      bloom, ht);
#endif

  if (opt_differences > 1)
    {
#ifdef SMART

#if 1
      generate_variants_1(hash,
                          sequence, seqlen,
                          v_gene, d_gene,
                          variant_list, variant_count,
                          bloom, ht);
#endif

      generate_variants_2_smart(hash,
                                sequence, seqlen,
                                v_gene, d_gene,
                                variant_list, variant_count,
                                bloom, ht);
#else
      generate_variants_2_all(hash,
                              sequence, seqlen,
                              v_gene, d_gene,
                              variant_list, variant_count,
                              bloom, ht);
#endif
    }
  else if (opt_differences > 0)
    {
      generate_variants_1(hash,
                          sequence, seqlen,
                          v_gene, d_gene,
                          variant_list, variant_count,
                          bloom, ht);
    }
  else if (opt_differences < 0)
    {
      /* dummy */
      generate_variants_2_smart(hash,
                                sequence, seqlen,
                                v_gene, d_gene,
                                variant_list, variant_count,
                                bloom, ht);
    }

#ifdef USEHASH
  hash_exit(ht);
#endif

#ifdef USEBLOOM
  bloom_exit(bloom);
#endif

#ifdef HASHSTATS
  printf("Hashstats: %" PRIu64 "\n", tries);
#endif

#ifdef DUMP
  exit(1);
#endif
}
