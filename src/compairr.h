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

#include <inttypes.h>

#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
#define PRIu64 "lu"
#endif
#endif

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
#define PRId64 "ld"
#endif
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <regex.h>
#include <limits.h>
#include <stdarg.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <chrono>
#include <string>
#include <vector>
#include <map>

#ifdef __APPLE__
#include <sys/resource.h>
#include <sys/sysctl.h>
#elif defined _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef __aarch64__

#include <arm_neon.h>

#elif defined __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

#define CAST_m128i_ptr(x) (reinterpret_cast<__m128i*>(x))

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#else

#error Unknown architecture
#endif

static_assert(INT_MAX > 32767, "Your compiler uses very short integers.");

/* constants */

#define PROG_CMD "compairr"
#define PROG_NAME "CompAIRR"
#define PROG_VERSION "1.3.0"
#define PROG_BRIEF "Comparison of Adaptive Immune Receptor Repertoires"

const unsigned int MAX_THREADS = 256;

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

extern int alphabet_size;

enum
  {
    summands_product,
    summands_ratio,
    summands_min,
    summands_max,
    summands_mean
  };

/* common data */

extern bool opt_alternative;
extern bool opt_cluster;
extern bool opt_help;
extern bool opt_ignore_counts;
extern bool opt_ignore_genes;
extern bool opt_indels;
extern bool opt_matrix;
extern bool opt_nucleotides;
extern bool opt_version;
extern char * opt_log;
extern char * opt_output_file;
extern char * opt_pairs;
extern char * opt_summands_string;
extern int64_t opt_differences;
extern int64_t opt_summands_int;
extern int64_t opt_threads;

extern FILE * outfile;
extern FILE * logfile;
extern FILE * pairsfile;

/* header files */

#include "util.h"
#include "arch.h"
#include "bloompat.h"
#include "cluster.h"
#include "db.h"
#include "hashtable.h"
#include "overlap.h"
#include "threads.h"
#include "variants.h"
#include "zobrist.h"
