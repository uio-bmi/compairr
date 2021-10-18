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
bool opt_cluster;
bool opt_existence;
bool opt_help;
bool opt_ignore_counts;
bool opt_ignore_genes;
bool opt_indels;
bool opt_matrix;
bool opt_nucleotides;
bool opt_version;
char * opt_log;
char * opt_output;
char * opt_pairs;
char * opt_summands_string;
int64_t opt_differences;
int64_t opt_summands_int;
int64_t opt_threads;

/* Other variables */

FILE * outfile = nullptr;
FILE * logfile = nullptr;
FILE * pairsfile = nullptr;

int alphabet_size;

static char dash[] = "-";
static char * DASH_FILENAME = dash;

static const char * summand_options[] = { "Product", "Ratio", "Min", "Max", "Mean" };

int64_t args_long(char * str, const char * option);
void args_show();
void args_usage();
void show_header();
void args_init(int argc, char **argv);
void open_files();
void close_files();

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
    fprintf(logfile, "Command (m/c/x):   Overlap (-m)\n");
  if (opt_cluster)
    fprintf(logfile, "Command (m/c/x):   Cluster (-c)\n");
  if (opt_existence)
    fprintf(logfile, "Command (m/c/x):   Existence (-x)\n");

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
  fprintf(logfile, "Threads (t):       %" PRId64 "\n", opt_threads);
  fprintf(logfile, "Output file (o):   %s\n", opt_output);
  if (opt_matrix || opt_existence)
    {
      fprintf(logfile, "Output format (a): %s\n", opt_alternative ? "Column" : "Matrix");
      fprintf(logfile, "Summands (s):      %s\n", summand_options[opt_summands_int]);
      fprintf(logfile, "Pairs file (p):    %s\n", opt_pairs ? opt_pairs : "(none)");
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
  fprintf(stderr, " -x, --existence             check existence of repertoire sequences in set\n");
  fprintf(stderr, " -c, --cluster               cluster sequences in one repertoire\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "General options:\n");
  fprintf(stderr, " -d, --differences INTEGER   number (0-2) of differences accepted (0)\n");
  fprintf(stderr, " -i, --indels                allow insertions or deletions\n");
  fprintf(stderr, " -f, --ignore-counts         ignore duplicate_count information\n");
  fprintf(stderr, " -g, --ignore-genes          ignore V and J gene information\n");
  fprintf(stderr, " -n, --nucleotides           compare nucleotides, not amino acids\n");
  fprintf(stderr, " -s, --summands STRING       sum product (default), ratio, min, max, or mean\n");
  fprintf(stderr, " -t, --threads INTEGER       number of threads to use (1)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output options:\n");
  fprintf(stderr, " -a, --alternative           output matrix results in column format\n");
  fprintf(stderr, " -p, --pairs FILENAME        output matching pairs to file (none)\n");
  fprintf(stderr, " -l, --log FILENAME          log to file (stderr)\n");
  fprintf(stderr, " -o, --output FILENAME       output results to file (stdout)\n");
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

  opt_help = false;
  opt_version = false;
  opt_matrix = false;
  opt_existence = false;
  opt_cluster = false;
  opt_differences = 0;
  opt_indels = false;
  opt_ignore_counts = false;
  opt_ignore_genes = false;
  opt_nucleotides = false;
  opt_summands_int = 0;
  opt_summands_string = NULL;
  opt_threads = 1;
  opt_alternative = false;
  opt_pairs = nullptr;
  opt_log = nullptr;
  opt_output = DASH_FILENAME;

  opterr = 1;

  char short_options[] = "acd:fghil:mno:p:s:t:vx";

  /* unused short option letters: bejkqruwxyz */

  static struct option long_options[] =
  {
    {"alternative",      no_argument,       nullptr, 'a' },
    {"cluster",          no_argument,       nullptr, 'c' },
    {"differences",      required_argument, nullptr, 'd' },
    {"ignore-counts",    no_argument,       nullptr, 'f' },
    {"ignore-genes",     no_argument,       nullptr, 'g' },
    {"help",             no_argument,       nullptr, 'h' },
    {"indels",           no_argument,       nullptr, 'i' },
    {"log",              required_argument, nullptr, 'l' },
    {"matrix",           no_argument,       nullptr, 'm' },
    {"nucleotides",      no_argument,       nullptr, 'n' },
    {"output",           required_argument, nullptr, 'o' },
    {"pairs",            required_argument, nullptr, 'p' },
    {"summands",         required_argument, nullptr, 's' },
    {"threads",          required_argument, nullptr, 't' },
    {"version",          no_argument,       nullptr, 'v' },
    {"existence",        no_argument,       nullptr, 'x' },
    {nullptr,            0,                 nullptr, 0   }
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
        /* summands */
        opt_summands_string = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'v':
        /* version */
        opt_version = true;
        break;

      case 'x':
        /* existence */
        opt_existence = true;
        break;

      default:
        show_header();
        args_usage();
        exit(1);
    }
  }

  int cmd_count = opt_help + opt_version + opt_matrix + opt_cluster + opt_existence;
  if (cmd_count == 0)
    fatal("Please specify a command (-h, -v, -c, -m or -x)");
  if (cmd_count > 1)
    fatal("Please specify just one command (-h, -v, -c, -m or -x)");

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
  else if (opt_cluster)
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

  if ((opt_threads < 1) || (opt_threads > MAX_THREADS))
    {
      fprintf(stderr, "\nError: Illegal number of threads specified with "
              "-t or --threads, must be in the range 1 to %u.\n", MAX_THREADS);
      exit(1);
    }

  if ((opt_differences < 0) || (opt_differences > 2))
    fatal("Differences specifed with -d or -differences must be 0, 1 or 2.");

  if (opt_indels && (opt_differences != 1))
    fatal("Indels are only allowed when d=1");

  if (opt_cluster)
    {
      if (opt_pairs)
        fatal("Option -p or --pairs is not allowed with -c or --cluster");
      if (opt_alternative)
        fatal("Option -a or --alternative is not allowed with -c or --cluster");
    }

  if (opt_summands_string)
    {
      opt_summands_int = -1;
      for(int i = 0; i < 5; i++)
        if (strcasecmp(opt_summands_string, summand_options[i]) == 0)
          {
            opt_summands_int = i;
            break;
          }
      if (opt_summands_int < 0)
        {
          fatal("Argument to -s or --summands must be product, ratio, min, max or mean");
        }
    }

  if (opt_nucleotides)
    alphabet_size = 4;
  else
    alphabet_size = 20;
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
  else
    cluster(input1_filename);

  show_time("End time:          ");

  close_files();
}
