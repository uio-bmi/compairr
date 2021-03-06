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

extern uint64_t * zobrist_tab_base;

inline uint64_t zobrist_value(unsigned int pos, unsigned char x)
{
  return zobrist_tab_base[alphabet_size * pos + x];
}

void zobrist_init(unsigned int longest,
                  unsigned int v_genes,
                  unsigned int d_genes);

void zobrist_exit();

uint64_t zobrist_hash(unsigned char * s,
                      unsigned int len,
                      int v_gene,
                      int d_gene);

uint64_t zobrist_hash_delete_first(unsigned char * s,
                                   unsigned int len,
                                   int v_gene,
                                   int d_gene);

uint64_t zobrist_hash_delete_first_two(unsigned char * s,
                                       unsigned int len,
                                       int v_gene,
                                       int d_gene);

uint64_t zobrist_hash_insert_first(unsigned char * s,
                                   unsigned int len,
                                   int v_gene,
                                   int d_gene);

uint64_t zobrist_hash_insert_first_two(unsigned char * s,
                                       unsigned int len,
                                       int v_gene,
                                       int d_gene);
