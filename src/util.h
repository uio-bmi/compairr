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

/* functions in util.cc */

int64_t gcd(int64_t a, int64_t b);
[[ noreturn ]] void fatal(const char * msg);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
void xfree(void * ptr);
char * xstrdup(const char * s);
void progress_init(const char * prompt, uint64_t size);
void progress_update(uint64_t progress);
void progress_done();
FILE * fopen_input(const char * filename);
FILE * fopen_output(const char * filename);
