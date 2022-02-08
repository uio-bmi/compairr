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

/* Variant information */

enum mutation_kind_enum
  {
   identical,
   substitution,
   deletion,
   insertion,
   sub_sub
  };

/*
  Info about the variants:

  For single mutations, only position1 and residue1 is used.

  For double mutations, position1 and residue1 applies to the first
  kind of mutation, while position2 and residue2 applies to the second
  kind of mutation. It is the order of the two mutations that matter,
  not their positions in the sequence.

  The residues residue1 and residue2 are the new residues that are
  present in the target sequence after insertion or substitution.
  For deletions, the residue is undefined.

  The positions position1 and position2 are the positions of the
  deletions, insertions or substitutions. The position1 value is
  the position of the first kind of mutation, while position2 is the
  position of the second kind of mutation. Note that
  position1 > position2 is not unusual.

  Positions start at zero.

  Positions of insertions refer to the new sequence.
  Positions of deletions refer to the original sequence.
  Positions of substitutions refer to the new sequence.

  For mutations of the same kind (del_del, sub_sub, ins_ins), the
  positions (position1 and position2) will always be in order;
  position1 < position2.

*/

struct var_s
{
  uint64_t hash;
  enum mutation_kind_enum kind;
  unsigned int pos1;
  unsigned int pos2;
  unsigned char residue1;
  unsigned char residue2;
};

void generate_variant_sequence(unsigned char * seed_sequence,
                               unsigned int seed_seqlen,
                               struct var_s * var,
                               unsigned char * seq,
                               unsigned int * seqlen);

bool check_variant(unsigned char * seed_sequence,
                   unsigned int seed_seqlen,
                   struct var_s * var,
                   unsigned char * amp_sequence,
                   unsigned int amp_seqlen);

void generate_variants(uint64_t hash,
                       unsigned char * sequence,
                       unsigned int seqlen,
                       uint64_t v_gene,
                       uint64_t d_gene,
                       struct var_s * variant_list,
                       unsigned int * variant_count);

uint64_t max_variants(uint64_t longest);
