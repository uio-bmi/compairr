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

/*

  This program uses Frederic Mahe's idea for swarm (d=1) to
  enumerate all variants of a sequence containing a single
  change (substitution, deletion or insertion) to quickly
  identify neighbour sequences using a hashing strategy.

  Please see the following publications for details:

  Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)
  Swarm: robust and fast clustering method for amplicon-based studies
  PeerJ 2:e593 https://doi.org/10.7717/peerj.593

  Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)
  Swarm v2: highly-scalable and high-resolution amplicon clustering
  PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420

*/

#include "compairr.h"

/* OPTIONS */

static char * progname;
static char * input1_filename;
static char * input2_filename;

bool opt_alternative;
bool opt_cdr3;
bool opt_cluster;
bool opt_distance;
bool opt_existence;
bool opt_help;
bool opt_ignore_counts;
bool opt_ignore_genes;
bool opt_ignore_unknown;
bool opt_indels;
bool opt_matrix;
bool opt_nucleotides;
bool opt_version;
bool opt_deduplicate;
char * opt_keep_columns;
char * opt_log;
char * opt_output;
char * opt_pairs;
char * opt_score_string;
int64_t opt_differences;
int64_t opt_score_int;
int64_t opt_threads;

/* Other variables */

const char * seq_header = nullptr;

FILE * outfile = nullptr;
FILE * logfile = nullptr;
FILE * pairsfile = nullptr;

int keep_columns_count = 0;
int * keep_columns_no = nullptr;
char ** keep_columns_names = nullptr;
char ** keep_columns_strings = nullptr;

int alphabet_size;

static char dash[] = "-";
static char * DASH_FILENAME = dash;

static const char * score_options[] =
  { "Product", "Ratio", "Min", "Max", "Mean", "MH", "Jaccard" };

static const char * score_descr[] =
  {
    "Sum of products of counts",
    "Sum of ratios of counts",
    "Sum of minimum of counts",
    "Sum of maximum of counts",
    "Sum of mean of counts",
    "Morisita-Horn index",
    "Jaccard index"
  };

int64_t args_long(char * str, const char * option);
void args_show();
void args_usage();
void show_header();
void args_init(int argc, char **argv);
void open_files();
void close_files();

bool parse_keep_columns()
{
  unsigned int len = strlen(opt_keep_columns);
  keep_columns_count = 1;
  for (unsigned int i = 0; i < len; i++)
    if (opt_keep_columns[i] == ',')
      keep_columns_count++;

  keep_columns_no = (int *) xmalloc
    (keep_columns_count * sizeof(int));

  keep_columns_names = (char **) xmalloc
    (keep_columns_count * sizeof(char *));

  keep_columns_strings = (char **) xmalloc
    (keep_columns_count * sizeof(char *));

  for (int j = 0; j < keep_columns_count; j++)
    keep_columns_no[j] = 0;

  keep_columns_count = 0;
  unsigned int curlen = 0;
  for (unsigned int i = 0; i < len; i++)
    {
      char c = opt_keep_columns[i];
      if (c == ',')
        {
          if (curlen == 0)
            return false;
          else
            {
              opt_keep_columns[i] = 0;
              keep_columns_names[keep_columns_count] =
                xstrdup(opt_keep_columns + i - curlen);
              opt_keep_columns[i] = ',';
              keep_columns_count++;
              curlen = 0;
            }
        }
      else if (((c >= 'A') && (c <= 'Z')) ||
               ((c >= 'a') && (c <= 'z')) ||
               ((c >= '0') && (c <= '9')) ||
               (c == '_'))
        {
          curlen++;
        }
      else
        {
          return false;
        }
    }

  if (curlen == 0)
    return false;

  keep_columns_names[keep_columns_count] =
    xstrdup(opt_keep_columns + len - curlen);
  keep_columns_count++;
  return true;
}

int64_t args_long(char * str, const char * option)
{
  char * endptr;
  int64_t temp = strtol(str, & endptr, 10);
  if (*endptr)
    {
      fprintf(stderr, "\nInvalid numeric argument for option %s\n", option);
      exit(1);
    }
  return temp;
}

void show_time(const char * prompt)
{
  const int time_string_max = 100;
  char time_string[time_string_max];
  const time_t clock = time(nullptr);
  const struct tm * timeptr = localtime(& clock);
  size_t time_string_len = strftime(time_string,
                                    time_string_max,
                                    "%a %b %d %T %Z %Y",
                                    timeptr);
  fprintf(logfile, "%s%s\n", prompt, time_string_len > 0 ? time_string : "?");
}

void args_show()
{
  if (opt_matrix)
    fprintf(logfile, "Command:           Overlap (-m)\n");
  if (opt_cluster)
    fprintf(logfile, "Command:           Cluster (-c)\n");
  if (opt_existence)
    fprintf(logfile, "Command:           Existence (-x)\n");
  if (opt_deduplicate)
    fprintf(logfile, "Command:           Deduplicate (--deduplicate)\n");

  if (opt_matrix)
    fprintf(logfile, "Repertoire set 1:  %s\n", input1_filename);
  else
    fprintf(logfile, "Repertoire:        %s\n", input1_filename);
  if (opt_matrix)
    fprintf(logfile, "Repertoire set 2:  %s\n", input2_filename ? input2_filename : "(same as set 1)");
  if (opt_existence)
    fprintf(logfile, "Repertoire set:    %s\n", input2_filename);

  fprintf(logfile, "Nucleotides (n):   %s\n", opt_nucleotides ? "Yes" : "No");
  fprintf(logfile, "Differences (d):   %" PRId64 "\n", opt_differences);
  fprintf(logfile, "Indels (i):        %s\n", opt_indels ? "Yes" : "No");
  fprintf(logfile, "Ignore counts (f): %s\n",
          opt_ignore_counts ? "Yes" : "No");
  fprintf(logfile, "Ignore genes (g):  %s\n",
          opt_ignore_genes ? "Yes" : "No");
  fprintf(logfile, "Ign. unknown (u):  %s\n",
          opt_ignore_unknown ? "Yes" : "No");
  fprintf(logfile, "Use cdr3 column:   %s\n",
          opt_cdr3 ? "Yes" : "No");
  fprintf(logfile, "Threads (t):       %" PRId64 "\n", opt_threads);
  fprintf(logfile, "Output file (o):   %s\n", opt_output);
  if (opt_matrix || opt_existence)
    {
      fprintf(logfile, "Output format (a): %s\n", opt_alternative ? "Column" : "Matrix");
      fprintf(logfile, "Score (s):         %s\n", score_descr[opt_score_int]);
      fprintf(logfile, "Pairs file (p):    %s\n", opt_pairs ? opt_pairs : "(none)");
      fprintf(logfile, "Keep columns:      %s\n", opt_keep_columns ? opt_keep_columns : "");
    }
  fprintf(logfile, "Log file (l):      %s\n", opt_log ? opt_log : "(stderr)");
}

void args_usage()
{
  fprintf(stderr, "Usage: %s [OPTIONS] TSVFILE1 [TSVFILE2]\n", PROG_CMD);
  fprintf(stderr, "\n");
  fprintf(stderr, "Commands:\n");
  fprintf(stderr, " -h, --help                  display this help and exit\n");
  fprintf(stderr, " -v, --version               display version information\n");
  fprintf(stderr, " -m, --matrix                compute overlap matrix between two sets\n");
  fprintf(stderr, " -x, --existence             check existence of sequences in repertoires\n");
  fprintf(stderr, " -c, --cluster               cluster sequences in one repertoire\n");
  fprintf(stderr, " -z, --deduplicate           deduplicate sequences in repertoires\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "General options:\n");
  fprintf(stderr, " -d, --differences INTEGER   number of differences accepted (0*)\n");
  fprintf(stderr, " -i, --indels                allow insertions or deletions when d=1\n");
  fprintf(stderr, " -f, --ignore-counts         ignore duplicate_count information\n");
  fprintf(stderr, " -g, --ignore-genes          ignore V and J gene information\n");
  fprintf(stderr, " -n, --nucleotides           compare nucleotides, not amino acids\n");
  fprintf(stderr, " -s, --score STRING          MH, Jaccard, product*, ratio, min, max, or mean\n");
  fprintf(stderr, " -t, --threads INTEGER       number of threads to use (1*-256)\n");
  fprintf(stderr, " -u, --ignore-unknown        ignore sequences with unknown symbols\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output options:\n");
  fprintf(stderr, " -a, --alternative           output results in three-column format, not matrix\n");
  fprintf(stderr, "     --cdr3                  use the cdr3(_aa) column instead of junction(_aa)\n");
  fprintf(stderr, "     --distance              include sequence distance in pairs file\n");
  fprintf(stderr, " -k, --keep-columns STRING   comma-separated columns to copy to pairs file\n");
  fprintf(stderr, " -l, --log FILENAME          log to file (stderr*)\n");
  fprintf(stderr, " -o, --output FILENAME       output results to file (stdout*)\n");
  fprintf(stderr, " -p, --pairs FILENAME        output matching pairs to file (none*)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "                             * default value\n");
  fprintf(stderr, "\n");
}

void show_header()
{
  fprintf(logfile, "%s %s - %s\n", PROG_NAME, PROG_VERSION, PROG_BRIEF);
  fprintf(logfile, "https://github.com/uio-bmi/compairr\n");
  fprintf(logfile, "\n");
}

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];
  input1_filename = nullptr;
  input2_filename = nullptr;

  opt_alternative = false;
  opt_cdr3 = false;
  opt_cluster = false;
  opt_deduplicate = false;
  opt_distance = false;
  opt_differences = 0;
  opt_existence = false;
  opt_help = false;
  opt_ignore_counts = false;
  opt_ignore_genes = false;
  opt_ignore_unknown = false;
  opt_indels = false;
  opt_keep_columns = nullptr;
  opt_log = nullptr;
  opt_matrix = false;
  opt_nucleotides = false;
  opt_output = DASH_FILENAME;
  opt_pairs = nullptr;
  opt_score_int = 0;
  opt_score_string = NULL;
  opt_threads = 1;
  opt_version = false;

  opterr = 1;

  char short_options[] = "acd:fghik:l:mno:p:s:t:uvxz";

  /* unused short option letters: bejqrwy */

  static struct option long_options[] =
  {
    {"alternative",      no_argument,       nullptr, 'a' },
    {"cdr3",             no_argument,       nullptr, 0   },
    {"cluster",          no_argument,       nullptr, 'c' },
    {"differences",      required_argument, nullptr, 'd' },
    {"distance",         no_argument,       nullptr, 0   },
    {"ignore-counts",    no_argument,       nullptr, 'f' },
    {"ignore-genes",     no_argument,       nullptr, 'g' },
    {"help",             no_argument,       nullptr, 'h' },
    {"indels",           no_argument,       nullptr, 'i' },
    {"keep-columns",     required_argument, nullptr, 'k' },
    {"log",              required_argument, nullptr, 'l' },
    {"matrix",           no_argument,       nullptr, 'm' },
    {"nucleotides",      no_argument,       nullptr, 'n' },
    {"output",           required_argument, nullptr, 'o' },
    {"pairs",            required_argument, nullptr, 'p' },
    {"score",            required_argument, nullptr, 's' },
    {"summands",         required_argument, nullptr, 's' },
    {"threads",          required_argument, nullptr, 't' },
    {"ignore-unknown",   no_argument,       nullptr, 'u' },
    {"version",          no_argument,       nullptr, 'v' },
    {"existence",        no_argument,       nullptr, 'x' },
    {"deduplicate",      no_argument,       nullptr, 'z' },
    {nullptr,            0,                 nullptr, 0   }
  };

  enum
    {
      option_alternative,
      option_cdr3,
      option_cluster,
      option_differences,
      option_distance,
      option_ignore_counts,
      option_ignore_genes,
      option_help,
      option_indels,
      option_keep_columns,
      option_log,
      option_matrix,
      option_nucleotides,
      option_output,
      option_pairs,
      option_score,
      option_summands,
      option_threads,
      option_ignore_unknown,
      option_version,
      option_existence,
      option_deduplicate
    };

  int used_options[26] = { 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0 };

  int option_index = 0;
  int c;

  while ((c = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1)
  {

    /* check if any option is specified more than once */

    if ((c >= 'a') && (c <= 'z'))
      {
        int optindex = c - 'a';
        if (used_options[optindex] == 1)
          {
            int longoptindex = 0;
            while (long_options[longoptindex].name)
              {
                if (long_options[longoptindex].val == c)
                  break;
                longoptindex++;
              }

            fprintf(stderr,
                    "Error: Option -%c or --%s specified more than once.\n",
                    c,
                    long_options[longoptindex].name);
            exit(1);
          }
        used_options[optindex] = 1;
      }

    switch(c)
      {
      case 'a':
        /* alternative */
        opt_alternative = true;
        break;

      case 'c':
        /* cluster */
        opt_cluster = true;
        break;

      case 'd':
        /* differences */
        opt_differences = args_long(optarg, "-d or --differences");
        break;

      case 'f':
        /* ignore-counts */
        opt_ignore_counts = true;
        break;

      case 'g':
        /* ignore-genes */
        opt_ignore_genes = true;
        break;

      case 'h':
        /* help */
        opt_help = true;
        break;

      case 'i':
        /* indels */
        opt_indels = true;
        break;

      case 'k':
        /* keep_columns */
        opt_keep_columns = optarg;
        break;

      case 'l':
        /* log */
        opt_log = optarg;
        break;

      case 'm':
        /* matrix */
        opt_matrix = true;
        break;

      case 'n':
        /* nucleotides */
        opt_nucleotides = true;
        break;

      case 'o':
        /* output-file */
        opt_output = optarg;
        break;

      case 'p':
        /* pairs-file */
        opt_pairs = optarg;
        break;

      case 's':
        /* score, summands */
        opt_score_string = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'u':
        /* ignore-unknown */
        opt_ignore_unknown = true;
        break;

      case 'v':
        /* version */
        opt_version = true;
        break;

      case 'x':
        /* existence */
        opt_existence = true;
        break;

      case 'z':
        /* deduplicate */
        opt_deduplicate = true;
        break;

      case 0:
        /* long options only */

        switch (option_index)
          {
          case option_cdr3:
            /* cdr3 */
            opt_cdr3 = true;
            break;

          case option_distance:
            /* distance */
            opt_distance = true;
            break;

          default:
            show_header();
            args_usage();
            exit(1);
          }
        break;

      default:
        show_header();
        args_usage();
        exit(1);
    }
  }

  int cmd_count = opt_help + opt_version + opt_matrix + opt_cluster + opt_existence + opt_deduplicate;
  if (cmd_count == 0)
    fatal("Please specify a command (--help, --version, --matrix, --existence, --cluster, or --deduplicate)");
  if (cmd_count > 1)
    fatal("Please specify just one command (--help, --version, --matrix, --existence, --cluster, or --deduplicate)");

  if (opt_help || opt_version)
    {
      if (optind != argc)
        fatal("Incorrect number of arguments");
    }
  else if (opt_matrix)
    {
      if (optind + 2 == argc)
        {
          input1_filename = argv[optind];
          input2_filename = argv[optind + 1];
        }
      else if (optind + 1 == argc)
        {
          input1_filename = argv[optind];
          input2_filename = 0;
        }
      else
        {
          fatal("Incorrect number of arguments. One or two input files must be specified.");
        }
    }
  else if (opt_existence)
    {
      if (optind + 2 == argc)
        {
          input1_filename = argv[optind];
          input2_filename = argv[optind + 1];
        }
      else
        {
          fatal("Incorrect number of arguments. Two input files must be specified.");
        }
    }
  else if (opt_cluster || opt_deduplicate)
    {
      if (optind + 1 == argc)
        {
          input1_filename = argv[optind];
        }
      else
        {
          fatal("Incorrect number of arguments. One input file must be specified.");
        }
    }

  if (opt_deduplicate)
    {
      if (opt_differences != 0)
        fatal("Option -d or --differences must be 0 for deduplication.");
      if (opt_indels)
        fatal("Option -i or --indels is not allowed for deduplication.");
    }

  if (opt_keep_columns)
    {
      if (! opt_pairs)
        fatal("Option --keep-columns only allowed with --pairs options.");
      if (! parse_keep_columns())
        fatal("Illegal list of columns with --keep-columns option. It must be a comma-separated list of column names. Allowed symbols: A-Z, a-z, _, and 0-9.");
    }

  if ((opt_threads < 1) || (opt_threads > MAX_THREADS))
    {
      fprintf(stderr, "\nError: Illegal number of threads specified with "
              "-t or --threads, must be in the range 1 to %u.\n", MAX_THREADS);
      exit(1);
    }

  if (opt_differences < 0)
    fatal("Differences specified with -d or -differences cannot be negative.");

  if (opt_indels && (opt_differences != 1))
    fatal("Indels are only allowed when d=1");

  if (opt_cluster)
    {
      if (opt_pairs)
        fatal("Option -p or --pairs is not allowed with -c or --cluster");
      if (opt_alternative)
        fatal("Option -a or --alternative is not allowed with -c or --cluster");
      if (opt_score_string)
        fatal("Option -s or --score is not allowed with -c or --cluster");
    }

  if (opt_score_string)
    {
      opt_score_int = -1;
      for(int i = 0; i < score_end; i++)
        if (strcasecmp(opt_score_string, score_options[i]) == 0)
          {
            opt_score_int = i;
            break;
          }
      if (opt_score_int < 0)
        {
          fatal("Argument to -s or --score must be MH, Jaccard, product, ratio, min, max or mean");
        }
    }

  if (! opt_matrix)
    {
      if (opt_score_int == score_mh)
        {
          fatal("The Morisita-Horn index is only allowed when computing repertoire overlap");
        }
      if (opt_score_int == score_jaccard)
        {
          fatal("The Jaccard index is only allowed when computing repertoire overlap");
        }
    }

  if (opt_differences > 0)
    {
      if (opt_score_int == score_mh)
        {
          fatal("The Morisita-Horn index is not defined when d>0");
        }
      if (opt_score_int == score_jaccard)
        {
          fatal("The Jaccard index is not defined when d>0");
        }
    }

  if (opt_nucleotides)
    alphabet_size = 4;
  else
    alphabet_size = 20;

  if (opt_cdr3)
    if (opt_nucleotides)
      seq_header = "cdr3";
    else
      seq_header = "cdr3_aa";
  else
    if (opt_nucleotides)
      seq_header = "junction";
    else
      seq_header = "junction_aa";
}

void open_files()
{
  /* open files */

  if (opt_log)
    {
      logfile = fopen_output(opt_log);
      if (! logfile)
        fatal("Unable to open log file for writing.");
    }

  outfile = fopen_output(opt_output);
  if (! outfile)
    fatal("Unable to open output file for writing.");

  if (opt_pairs)
    {
      pairsfile = fopen_output(opt_pairs);
      if (! pairsfile)
        fatal("Unable to open pairs file for writing.");
    }
}

void close_files()
{
  if (pairsfile)
    fclose(pairsfile);

  if (outfile)
    fclose(outfile);

  if (logfile)
    fclose(logfile);
}

int main(int argc, char** argv)
{
  logfile = stderr;

  arch_srandom(1);

  args_init(argc, argv);

  open_files();

  if (opt_version || opt_help)
    {
      show_header();
      if (opt_help)
        args_usage();
      close_files();
      exit(0);
    }

  show_header();

  show_time("Start time:        ");

  args_show();

  fprintf(logfile, "\n");

  if (opt_matrix || opt_existence)
    overlap(input1_filename, input2_filename);
  else if (opt_deduplicate)
    dedup(input1_filename);
  else
    cluster(input1_filename);

  show_time("End time:          ");

  if (keep_columns_no)
    {
      xfree(keep_columns_no);
      keep_columns_no = nullptr;
    }

  if (keep_columns_names)
    {
      xfree(keep_columns_names);
      keep_columns_names = nullptr;
    }

  if (keep_columns_strings)
    {
      xfree(keep_columns_strings);
      keep_columns_strings = nullptr;
    }

  close_files();
}
