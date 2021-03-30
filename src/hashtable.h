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

struct hashtable_s
{
  uint64_t hash_mask;
  unsigned char * hash_occupied;
  uint64_t * hash_values;
  unsigned int * hash_data;
  uint64_t hash_tablesize;
};

inline uint64_t hash_get_tablesize(struct hashtable_s * ht)
{
  return ht->hash_tablesize;
}

inline uint64_t hash_getindex(struct hashtable_s * ht, uint64_t hash)
{
  // Shift bits right to get independence from the simple Bloom filter hash
  hash = hash >> 32;
  return hash & ht->hash_mask;
}

inline uint64_t hash_getnextindex(struct hashtable_s * ht, uint64_t j)
{
  return (j+1) & ht->hash_mask;
}

inline void hash_set_occupied(struct hashtable_s * ht, uint64_t j)
{
  ht->hash_occupied[j >> 3] |= (1 << (j & 7));
}

inline bool hash_is_occupied(struct hashtable_s * ht, uint64_t j)
{
  return ht->hash_occupied[j >> 3] & (1 << (j & 7));
}

inline void hash_set_value(struct hashtable_s * ht, uint64_t j, uint64_t hash)
{
  ht->hash_values[j] = hash;
}

inline bool hash_compare_value(struct hashtable_s * ht,
                               uint64_t j, uint64_t hash)
{
  return (ht->hash_values[j] == hash);
}

inline unsigned int hash_get_data(struct hashtable_s * ht, uint64_t j)
{
  return ht->hash_data[j];
}

inline void hash_set_data(struct hashtable_s * ht, uint64_t j, unsigned int x)
{
  ht->hash_data[j] = x;
}

void hash_zap(struct hashtable_s * ht);

struct hashtable_s * hash_init(uint64_t sequences);

void hash_exit(struct hashtable_s * ht);
