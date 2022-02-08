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

    d=2
    substitutions: (ab_size - 1) * (ab_size - 1) * longest * (longest - 1) / 2
  */

  uint64_t maxvar = 0;

  // d = 0
  // identical non-variant
  maxvar += 1;

  if (opt_differences >= 1)
    {
      // d = 1
      // substitutions
      maxvar += longest * (alphabet_size - 1);

      if (opt_indels)
        {
          // deletions
          maxvar += longest;

          // insertions
          maxvar += (longest + 1) * (alphabet_size - 1) + 1;
        }
    }

  if (opt_differences >= 2)
    {
      // d = 2
      // substitutions
      maxvar += longest * (longest - 1) / 2 *
        (alphabet_size - 1) * (alphabet_size - 1);

      if (opt_indels)
        {
          // deletions & insertions
          fatal("Indels not supported for d>1");
        }
    }

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
  /* make sure seed with given variant is really identical to amp */
  /* we know the hashes are identical */

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
                        unsigned char residue2)
{
  var_s * v = variant_list + (*variant_count)++;
  v->hash = hash;
  v->kind = kind;
  v->pos1 = pos1;
  v->residue1 = residue1;
  v->pos2 = pos2;
  v->residue2 = residue2;
}

void generate_variants_0(uint64_t hash,
                         var_s * variant_list,
                         unsigned int * variant_count)
{
  /* identical non-variant */
  add_variant(hash,
              variant_list, variant_count,
              identical, 0, 0, 0, 0);
}

void generate_variants_1(uint64_t hash,
                         unsigned char * sequence,
                         unsigned int seqlen,
                         uint64_t v_gene,
                         uint64_t d_gene,
                         var_s * variant_list,
                         unsigned int * variant_count)
{
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
                        substitution, i, v, 0, 0);
          }
    }

  /* indels */

  if (opt_indels)
    {
      /* deletions */

      if (seqlen > 1)
        {
          hash = zobrist_hash_delete_first
            (reinterpret_cast<unsigned char *> (sequence),
             seqlen,
             v_gene,
             d_gene);
          add_variant(hash,
                      variant_list, variant_count,
                      deletion, 0, 0, 0, 0);
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
                              deletion, i, 0, 0, 0);
                  deleted = v;
                }
            }
        }

      /* insertions */

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
                      insertion, 0, v, 0, 0);
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
                            insertion, i + 1, v, 0, 0);
              }
        }
    }
}

void generate_variants_2(uint64_t hash,
                         unsigned char * sequence,
                         unsigned int seqlen,
                         uint64_t v_gene,
                         uint64_t d_gene,
                         var_s * variant_list,
                         unsigned int * variant_count)
{
  (void) v_gene;
  (void) d_gene;

  /* generate all double substitutions */

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
                                      sub_sub, i, v, j, w);
                        }
                    }
                }
            }
        }
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
  generate_variants_0(hash,
                      variant_list, variant_count);

  if (opt_differences >= 1)
    {
      generate_variants_1(hash,
                          sequence, seqlen,
                          v_gene, d_gene,
                          variant_list, variant_count);
    }

  if (opt_differences >= 2)
    {
      generate_variants_2(hash,
                          sequence, seqlen,
                          v_gene, d_gene,
                          variant_list, variant_count);
    }
}
