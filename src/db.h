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

/* structures and data types */

struct db;



/* functions in db.cc */

void db_init();

void db_exit();

struct db * db_create();

void db_free(struct db * d);

void db_read(struct db * d,
             const char * filename,
             bool require_sequence_id,
             const char * default_repertoire_id);

uint64_t db_getsequencecount(struct db * d);

uint64_t db_get_repertoire_count(struct db * d);

uint64_t db_getresiduescount(struct db * d);

unsigned int db_getlongestsequence(struct db * d);

char * db_getsequence(struct db * d, uint64_t seqno);

unsigned int db_getsequencelen(struct db * d, uint64_t seqno);

uint64_t db_gethash(struct db * d, uint64_t seqno);

uint64_t db_get_v_gene(struct db * d, uint64_t seqno);

uint64_t db_get_j_gene(struct db * d, uint64_t seqno);

uint64_t db_get_count(struct db * d, uint64_t seqno);

int db_get_repertoire_id_no(struct db * d, uint64_t seqno);

const char * db_get_repertoire_id(struct db * d, int repertoire_id_no);

char * db_get_sequence_id(struct db * d, uint64_t seqno);

void db_hash(struct db * d);

uint64_t db_get_v_gene_count();

uint64_t db_get_j_gene_count();

const char * db_get_v_gene_name(struct db * d, uint64_t seqno);

const char * db_get_j_gene_name(struct db * d, uint64_t seqno);

void db_fprint_sequence(FILE * f, struct db * d, uint64_t seqno);

char * db_get_keep_columns(struct db * d, uint64_t seqno);
