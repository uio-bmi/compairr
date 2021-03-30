/*
    Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe

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

#include "vdjsearch.h"

/* OPTIONS */

static char * progname;
static char * input1_filename;
static char * input2_filename;

char * opt_log;
char * opt_output_file;

int64_t opt_differences;
int64_t opt_help;
int64_t opt_indels;
int64_t opt_threads;
int64_t opt_version;

/* Other variables */

FILE * outfile = nullptr;
FILE * logfile = nullptr;
FILE * fp_seeds = nullptr;
FILE * network_file = nullptr;

static char dash[] = "-";
static char * DASH_FILENAME = dash;

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

void args_show()
{
  fprintf(logfile, "Repertoire set 1:  %s\n", input1_filename);
  fprintf(logfile, "Repertoire set 2:  %s\n", input2_filename);
  fprintf(logfile, "Differences (d):   %lld\n", opt_differences);
  fprintf(logfile, "Indels (i):        %s\n", opt_indels ? "Yes" : "No");
  fprintf(logfile, "Output file (o):   %s\n", opt_output_file);
  fprintf(logfile, "Threads: (t)       %" PRId64 "\n", opt_threads);
}

void args_usage()
{
  fprintf(stderr, "Usage: %s [OPTIONS] TSVFILE1 TSVFILE2\n", PROG_NAME);
  fprintf(stderr, "\n");
  fprintf(stderr, "General options:\n");
  fprintf(stderr, " -d, --differences INTEGER           number (0-2) of differences accepted (1)\n");
  fprintf(stderr, " -i, --indels                        allow insertions or deletions (no)\n");
  fprintf(stderr, " -h, --help                          display this help and exit\n");
  fprintf(stderr, " -t, --threads INTEGER               number of threads to use (1)\n");
  fprintf(stderr, " -v, --version                       display version information and exit\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output options:\n");
  fprintf(stderr, " -l, --log FILENAME                  log to file, not to stderr\n");
  fprintf(stderr, " -o, --output-file FILENAME          output results to file (stdout)\n");
  fprintf(stderr, "\n");
}

void show_header()
{
  fprintf(logfile,
          "%s %s - %s\n\n",
          PROG_NAME,
          PROG_VERSION,
          PROG_BRIEF);
}

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];
  input1_filename = 0;
  input2_filename = 0;

  opt_differences = 1;
  opt_help = 0;
  opt_indels = 0;
  opt_log = nullptr;
  opt_output_file = DASH_FILENAME;
  opt_threads = 1;
  opt_version = 0;

  opterr = 1;

  char short_options[] = "d:hil:o:t:v";

  /* unused short option letters: abcefgjkmnpqrsuwxyz */

  static struct option long_options[] =
  {
    {"differences",           required_argument, nullptr, 'd' },
    {"help",                  no_argument,       nullptr, 'h' },
    {"indels",                no_argument,       nullptr, 'i' },
    {"log",                   required_argument, nullptr, 'l' },
    {"output-file",           required_argument, nullptr, 'o' },
    {"threads",               required_argument, nullptr, 't' },
    {"version",               no_argument,       nullptr, 'v' },
    {nullptr,                 0,                 nullptr, 0 }
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
      case 'd':
        /* differences */
        opt_differences = args_long(optarg, "-d or --differences");
        break;

      case 'h':
        /* help */
        opt_help = 1;
        break;

      case 'i':
        /* indels */
        opt_indels = 1;
        break;

      case 'l':
        /* log */
        opt_log = optarg;
        break;

      case 'o':
        /* output-file */
        opt_output_file = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'v':
        /* version */
        opt_version = 1;
        break;

      default:
        show_header();
        args_usage();
        exit(1);
    }
  }

  if (opt_help || opt_version)
    {
      if (optind != argc)
        fatal("Incorrect number of arguments");
    }
  else
    {
      if (optind + 2 == argc)
        {    
          input1_filename = argv[optind];
          input2_filename = argv[optind + 1];
        }
      else
        {
          fatal("Incorrect number of arguments (two input files must be specified)");
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

#if 0
  if (opt_indels && (opt_differences != 1))
    fatal("Indels is only allowed when d=1");
#endif
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

  outfile = fopen_output(opt_output_file);
  if (! outfile)
    fatal("Unable to open output file for writing.");

}

void close_files()
{
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

  args_show();

  fprintf(logfile, "\n");

  overlap(input1_filename, input2_filename);
  
  close_files();
}
