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

#include <string>
#include <map>
#include <vector>

/* How much memory for residues and sequences should we allocate each time? */

#define MEMCHUNK 1048576
#define SEQCHUNK 65536

static signed char map_aa[256] =
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1,  2,  3,  4,  5,  6,  7, -1,  8,  9, 10, 11, -1,
    12, 13, 14, 15, 16, -1, 17, 18, -1, 19, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1,  2,  3,  4,  5,  6,  7, -1,  8,  9, 10, 11, -1,
    12, 13, 14, 15, 16, -1, 17, 18, -1, 19, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

static signed char map_nt[256] =
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

const static char * aa_chars = "ACDEFGHIKLMNPQRSTVWY";
const static char * nt_chars = "acgt";
const char * EMPTYSTRING = "";

struct seqinfo_s
{
  uint64_t hash;
  uint64_t count;
  char * sequence_id;
  char * seq;
  char * keep; /* extra columns to keep, tab separated */
  unsigned int seqlen;
  int repertoire_id_no;
  int v_gene_no;
  int j_gene_no;
};

typedef struct seqinfo_s seqinfo_t;

struct db
{
  seqinfo_t * seqindex;
  uint64_t seqindex_alloc;
  uint64_t sequences;
  unsigned int longest;
  unsigned int shortest;
  char * residues_p;
  uint64_t residues_alloc;
  uint64_t residues_count;
  uint64_t total_duplicate_count;
  uint64_t repertoire_count;
  uint64_t ignored_unknown;
  uint64_t ignored_empty;
  std::vector<std::string> repertoire_id_vector;
  std::map<std::string, int> repertoire_id_map;
  int col_junction;
  int col_junction_aa;
  int col_cdr3;
  int col_cdr3_aa;
  int col_duplicate_count;
  int col_v_call;
  int col_j_call;
  int col_repertoire_id;
  int col_sequence_id;
};

/* v and j genes are common to both */

static std::vector<std::string> v_gene_vector;
static std::map<std::string, int> v_gene_map;

static std::vector<std::string> j_gene_vector;
static std::map<std::string, int> j_gene_map;

void db_init()
{
  v_gene_vector.clear();
  v_gene_map.clear();
  j_gene_vector.clear();
  j_gene_map.clear();
}

void db_exit()
{
  v_gene_vector.clear();
  v_gene_map.clear();
  j_gene_vector.clear();
  j_gene_map.clear();
}

struct db * db_create()
{
  struct db * d = new db;

  d->seqindex = nullptr;
  d->seqindex_alloc = 0;
  d->sequences = 0;
  d->longest = 0;
  d->shortest = UINT_MAX;
  d->residues_p = nullptr;
  d->residues_alloc = 0;
  d->residues_count = 0;
  d->total_duplicate_count = 0;
  d->repertoire_count = 0;
  d->repertoire_id_vector.clear();
  d->repertoire_id_map.clear();
  d->col_junction = 0;
  d->col_junction_aa = 0;
  d->col_cdr3 = 0;
  d->col_cdr3_aa = 0;
  d->col_duplicate_count = 0;
  d->col_v_call = 0;
  d->col_j_call = 0;
  d->col_repertoire_id = 0;
  d->col_sequence_id = 0;

  return d;
}

void parse_airr_tsv_header(char * line,
                           struct db * d,
                           bool require_sequence_id)
{
  char delim[] = "\t";
  char * string = line;
  char * token = nullptr;

  int i = 1;

  while ((token = strsep(& string, delim)) != nullptr)
    {
      if (strcmp(token, "repertoire_id") == 0)
        {
          d->col_repertoire_id = i;
        }
      else if (strcmp(token, "sequence_id") == 0)
        {
          d->col_sequence_id = i;
        }
      else if (strcmp(token, "duplicate_count") == 0)
        {
          d->col_duplicate_count = i;
        }
      else if (strcmp(token, "v_call") == 0)
        {
          d->col_v_call = i;
        }
      else if (strcmp(token, "j_call") == 0)
        {
          d->col_j_call = i;
        }
      else if (strcmp(token, "junction") == 0)
        {
          d->col_junction = i;
        }
      else if (strcmp(token, "junction_aa") == 0)
        {
          d->col_junction_aa = i;
        }
      else if (strcmp(token, "cdr3") == 0)
        {
          d->col_cdr3 = i;
        }
      else if (strcmp(token, "cdr3_aa") == 0)
        {
          d->col_cdr3_aa = i;
        }

      for (int j = 0; j < keep_columns_count; j++)
        {
          if (strcmp(token, keep_columns_names[j]) == 0)
            keep_columns_no[j] = i;
        }
      i++;
    }

  if (! (d->col_sequence_id     || ! require_sequence_id)   ||
      ! (d->col_duplicate_count ||   opt_ignore_counts)     ||
      ! (d->col_v_call          ||   opt_ignore_genes)      ||
      ! (d->col_j_call          ||   opt_ignore_genes)      ||
      ! (d->col_junction        || ! opt_nucleotides        || opt_cdr3   ) ||
      ! (d->col_junction_aa     ||   opt_nucleotides        || opt_cdr3   ) ||
      ! (d->col_cdr3            || ! opt_nucleotides        || ! opt_cdr3 ) ||
      ! (d->col_cdr3_aa         ||   opt_nucleotides        || ! opt_cdr3 ))
    {
      fprintf(logfile,
        "\nMissing essential column(s) in header of AIRR TSV input file:");

      if (require_sequence_id && (! d->col_sequence_id))
        fprintf(logfile, " sequence_id");
      if ((! opt_ignore_counts) && (! d->col_duplicate_count))
        fprintf(logfile, " duplicate_count");
      if (! opt_ignore_genes)
        {
          if (! d->col_v_call)
            fprintf(logfile, " v_call");
          if (! d->col_j_call)
            fprintf(logfile, " j_call");
        }
      if (opt_cdr3)
        {
          if (opt_nucleotides)
            {
              if (! d->col_cdr3)
                fprintf(logfile, " cdr3");
            }
          else
            {
              if (! d->col_cdr3_aa)
                fprintf(logfile, " cdr3_aa");
            }
        }
      else
        {
          if (opt_nucleotides)
            {
              if (! d->col_junction)
                fprintf(logfile, " junction");
            }
          else
            {
              if (! d->col_junction_aa)
                fprintf(logfile, " junction_aa");
            }
        }

      fprintf(logfile, "\n");
      exit(1);
    }

  bool any_missing = false;
  for (int j = 0; j < keep_columns_count; j++)
    if (keep_columns_no[j] < 1)
      any_missing = true;
  if (any_missing)
    {
      fprintf(logfile,
              "\nWarning: missing column(s) to keep in header:");
      for (int j = 0; j < keep_columns_count; j++)
        if (keep_columns_no[j] < 1)
          fprintf(logfile, " %s", keep_columns_names[j]);
      fprintf(logfile, "\n");
    }
}

void parse_airr_tsv_line(char * line,
                         uint64_t lineno,
                         struct db * d,
                         bool require_sequence_id,
                         const char * default_repertoire_id)
{
  const char * repertoire_id = nullptr;
  const char * sequence_id = nullptr;
  const char * duplicate_count = nullptr;
  const char * v_call = nullptr;
  const char * j_call = nullptr;
  const char * junction = nullptr;
  const char * junction_aa = nullptr;
  const char * cdr3 = nullptr;
  const char * cdr3_aa = nullptr;

  for (int k = 0; k < keep_columns_count; k++)
    keep_columns_strings[k] = nullptr;

  char delim[] = "\t";
  char * string = line;
  char * token = nullptr;

  int i = 1;

  while ((token = strsep(& string, delim)) != nullptr)
    {
      if (i == d->col_repertoire_id)
        {
          repertoire_id = token;
        }
      else if (i == d->col_sequence_id)
        {
          sequence_id = token;
        }
      else if (i == d->col_duplicate_count)
        {
          duplicate_count = token;
        }
      else if (i == d->col_v_call)
        {
          v_call = token;
        }
      else if (i == d->col_j_call)
        {
          j_call = token;
        }
      else if (i == d->col_junction)
        {
          junction = token;
        }
      else if (i == d->col_junction_aa)
        {
          junction_aa = token;
        }
      else if (i == d->col_cdr3)
        {
          cdr3 = token;
        }
      else if (i == d->col_cdr3_aa)
        {
          cdr3_aa = token;
        }

      for (int k = 0; k < keep_columns_count; k++)
        if (i == keep_columns_no[k])
          keep_columns_strings[k] = token;

      i++;
    }


  /* make room for another entry */

  if (d->sequences >= d->seqindex_alloc)
    {
      d->seqindex_alloc += SEQCHUNK;
      d->seqindex = static_cast<seqinfo_t *>
        (xrealloc(d->seqindex, d->seqindex_alloc * sizeof(seqinfo_s)));
    }

  seqinfo_t * p = d->seqindex + d->sequences;


  /* make room for more residues */

  unsigned int len_estimate = 0;
  if (opt_cdr3)
    {
      if (opt_nucleotides)
        len_estimate = strlen(cdr3);
      else
        len_estimate = strlen(cdr3_aa);
    }
  else
    {
      if (opt_nucleotides)
        len_estimate = strlen(junction);
      else
        len_estimate = strlen(junction_aa);
    }

  if (d->residues_count + len_estimate > d->residues_alloc)
    {
      d->residues_alloc += MEMCHUNK;
      d->residues_p = static_cast<char *>
        (xrealloc(d->residues_p, d->residues_alloc));
    }


  /* scan and store sequence */

  char * q = d->residues_p + d->residues_count;
  unsigned int seqlen = 0;
  bool ignore_seq = false;

  for(unsigned int i = 0; i < len_estimate; i++)
    {
      unsigned char c;
      signed char m;
      if (opt_nucleotides)
        {
          if (opt_cdr3)
            c = cdr3[i];
          else
            c = junction[i];
          m = map_nt[static_cast<unsigned int>(c)];
        }
      else
        {
          if (opt_cdr3)
            c = cdr3_aa[i];
          else
            c = junction_aa[i];
          m = map_aa[static_cast<unsigned int>(c)];
        }

      if (m >= 0)
        {
          *q++ = m;
          seqlen++;
        }
      else
        {
          if ((c >= 32) && (c <= 126))
            {
              if (opt_ignore_unknown)
                {
                  ignore_seq = true;
                  d->ignored_unknown++;
                }
              else
                {
                  fprintf(logfile,
                          "\n\nError: Illegal character '%c' in sequence "
                          "on line %" PRIu64 ". Use -u to ignore.\n",
                          c,
                          lineno);
                  exit(1);
                }
            }
          else
            {
              fprintf(logfile,
                      "\n\nError: Illegal character (ascii no %d) in sequence "
                      "on line %" PRIu64 "\n",
                      c,
                      lineno);
              exit(1);
            }
        }
    }

  if (seqlen == 0)
    {
      if (opt_ignore_empty)
        {
          ignore_seq = true;
          d->ignored_empty++;
        }
      else
        {
          fprintf(logfile,
                  "\n\nError: Empty sequence in sequence "
                  "on line %" PRIu64 ". Use -e to ignore.\n",
                  lineno);
          exit(1);
        }
    }

  if (ignore_seq)
    {
      return;
    }
  else
    {
      d->residues_count += seqlen;
      p->seqlen = seqlen;
      if (seqlen > d->longest)
        d->longest = seqlen;
      if (seqlen < d->shortest)
        d->shortest = seqlen;
    }


  /* handle repertoire_id */

  if (! repertoire_id)
    {
      repertoire_id = default_repertoire_id;
    }

  auto r_it = d->repertoire_id_map.find(repertoire_id);
  if (r_it != d->repertoire_id_map.end())
    {
      p->repertoire_id_no = r_it->second;
    }
  else
    {
      p->repertoire_id_no = d->repertoire_id_vector.size();
      d->repertoire_id_vector.push_back(repertoire_id);
      d->repertoire_id_map.insert({repertoire_id, p->repertoire_id_no});
    }


  /* handle sequence_id */

  if (sequence_id && *sequence_id)
    {
      p->sequence_id = xstrdup(sequence_id);
    }
  else if (require_sequence_id)
    {
      fprintf(logfile,
              "\n\nError: missing or empty sequence_id value on line %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }
  else
    {
      p->sequence_id = nullptr;
    }


  /* handle duplicate_count */

  if (duplicate_count && *duplicate_count)
    {
      char * endptr = nullptr;
      long count = strtol(duplicate_count, &endptr, 10);
      if (endptr && (*endptr == 0) && (count >= 1))
        {
          p->count = count;
        }
      else
        {
          fprintf(logfile, "\n\nError: Illegal duplicate_count on line %"
                  PRIu64 ": %s\n", lineno, duplicate_count);
          exit(1);
        }
    }
  else if (opt_ignore_counts)
    {
      p->count = 1;
    }
  else
    {
      fprintf(logfile,
              "\n\nError: missing or empty duplicate_count on line %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }

  d->total_duplicate_count += p->count;


  /* handle v_call */

  if (! opt_ignore_genes && ! (v_call && *v_call))
    {
      fprintf(logfile,
              "\n\nError: missing or empty v_call value on line %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }

  if (! v_call)
    {
      v_call = EMPTYSTRING;
    }

  auto v_it = v_gene_map.find(v_call);
  if (v_it != v_gene_map.end())
    {
      p->v_gene_no = v_it->second;
    }
  else
    {
      p->v_gene_no = v_gene_vector.size();
      v_gene_vector.push_back(v_call);
      v_gene_map.insert({v_call, p->v_gene_no});
    }


  /* handle j_call */

  if (! opt_ignore_genes && ! (j_call && *j_call))
    {
      fprintf(logfile,
              "\n\nError: missing or empty j_call value on line %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }

  if (! j_call)
    {
      j_call = EMPTYSTRING;
    }

  auto j_it = j_gene_map.find(j_call);
  if (j_it != j_gene_map.end())
    {
      p->j_gene_no = j_it->second;
    }
  else
    {
      p->j_gene_no = j_gene_vector.size();
      j_gene_vector.push_back(j_call);
      j_gene_map.insert({j_call, p->j_gene_no});
    }


  /* handle junction(_aa) or cdr3(_aa) */

  bool seq_ok = false;
  if (opt_nucleotides)
    {
      if (opt_cdr3)
        {
          seq_ok = cdr3 && *cdr3;
        }
      else
        {
          seq_ok = junction && *junction;
        }
    }
  else
    {
      if (opt_cdr3)
        {
          seq_ok = cdr3_aa && *cdr3_aa;
        }
      else
        {
          seq_ok = junction_aa && *junction_aa;
        }
    }

  if (! seq_ok)
    {
      fprintf(logfile,
              "\n\nError: missing or empty %s value on line %"
              PRIu64 "\n",
              seq_header,
              lineno);
      exit(1);
    }


  /* handle keep_columns */

  unsigned int len = 0;
  for (int k = 0; k < keep_columns_count; k++)
    {
      if (keep_columns_strings[k])
        len += strlen(keep_columns_strings[k]);
      len++;
    }
  if (len > 0)
    p->keep = (char *) xmalloc(len);
  else
    p->keep = nullptr;

  len = 0;
  bool first = true;
  for (int k = 0; k < keep_columns_count; k++)
    {
      if (first)
        first = false;
      else
        p->keep[len++] = '\t';
      if (keep_columns_strings[k])
        {
          strcpy(p->keep + len, keep_columns_strings[k]);
          len += strlen(keep_columns_strings[k]);
          keep_columns_strings[k] = nullptr;
        }
    }
  if (p->keep)
    p->keep[len] = 0;

  p->hash = 0;

  d->sequences++;
}

void db_read(struct db * d,
             const char * filename,
             bool require_sequence_id,
             const char * default_repertoire_id)
{
  FILE * fp = nullptr;
  if (filename)
    {
      fp = fopen_input(filename);
      if (!fp)
        {
          fprintf(logfile,
                  "\nError: Unable to open input data file (%s).\n",
                  filename);
          exit(1);
        }
    }
  else
    fp = stdin;

  /* get file size */

  struct stat fs;

  if (fstat(fileno(fp), & fs))
    {
      fprintf(logfile, "\nUnable to fstat on input file (%s)\n", filename);
      exit(1);
    }
  bool is_regular = S_ISREG(fs.st_mode);
  uint64_t filesize = is_regular ? (uint64_t)(fs.st_size) : 0;
  uint64_t fileread = 0;

  if (! is_regular)
    fprintf(logfile, "Waiting for data from standard input...\n");

  size_t line_alloc = 4096;
  char * line = (char *) xmalloc(line_alloc);
  uint64_t lineno = 0;
  ssize_t linelen = 0;

  d->longest = 0;
  d->shortest = UINT_MAX;
  d->ignored_unknown = 0;
  d->ignored_empty = 0;

  int state = 0;

  progress_init("Reading sequences:", filesize);

  linelen = getline(& line, & line_alloc, fp);

  if (linelen < 0)
    fatal("Unable to read from the input file");

  fileread += linelen;

  if ((linelen > 0) && (line[linelen-1] == '\n'))
    {
      line[linelen-1] = 0;
      linelen--;
    }

  if ((linelen > 0) && (line[linelen-1] == '\r'))
    {
      line[linelen-1] = 0;
      linelen--;
    }

  while (linelen >= 0)
    {
      lineno++;

      if (state == 0)
        {
          if (line[0] == '#')
            {
              /* ignore initial comment section */
            }
          else if (line[0] == '@')
            {
              /* ignore initial comment section */
            }
          else
            {
              parse_airr_tsv_header(line,
                                    d,
                                    require_sequence_id);
              state = 1;
            }
        }
      else
        {
          parse_airr_tsv_line(line,
                              lineno,
                              d,
                              require_sequence_id,
                              default_repertoire_id);
        }

      /* update progress */

      if (is_regular)
        progress_update(fileread);

      /* get next line */

      linelen = getline(& line, & line_alloc, fp);

      if (linelen < 0)
        break;

      fileread += linelen;

      /* remove LF at end of line */

      if ((linelen > 0) && (line[linelen-1] == '\n'))
        {
          line[linelen-1] = 0;
          linelen--;
        }

      /* remove CR at end of line if from DOS/Windows */

      if ((linelen > 0) && (line[linelen-1] == '\r'))
        {
          line[linelen-1] = 0;
          linelen--;
        }
    }

  progress_done();

  if (line)
    xfree(line);
  line = nullptr;

  fclose(fp);

  d->repertoire_count = d->repertoire_id_vector.size();

  if (d->ignored_unknown > 0)
    fprintf(logfile, "%" PRIu64 " sequences with unknown symbols ignored.\n", d->ignored_unknown);

  if (d->ignored_empty > 0)
    fprintf(logfile, "%" PRIu64 " empty sequences ignored.\n", d->ignored_empty);

  if (d->sequences > 0)
    {
      fprintf(logfile,
              "Repertoires:       %" PRIu64 "\n"
              "Sequences:         %" PRIu64 "\n"
              "Residues:          %" PRIu64 "\n"
              "Shortest:          %u\n"
              "Longest:           %u\n"
              "Average length:    %.1lf\n"
              "Total dupl. count: %" PRIu64 "\n",
              d->repertoire_count,
              d->sequences,
              d->residues_count,
              d->shortest,
              d->longest,
              1.0 * d->residues_count / d->sequences,
              d->total_duplicate_count);
    }
  else
    {
      fprintf(logfile,
              "Repertoires:       %" PRIu64 "\n"
              "Sequences:         %" PRIu64 "\n"
              "Residues:          %" PRIu64 "\n"
              "Shortest:          -\n"
              "Longest:           -\n"
              "Average length:    -\n"
              "Total dupl. count: %" PRIu64 "\n",
              d->repertoire_count,
              d->sequences,
              d->residues_count,
              d->total_duplicate_count);
    }

  /* add sequence pointers to index table */

  progress_init("Indexing:         ", d->sequences);
  char * r = d->residues_p;
  for(uint64_t i = 0; i < d->sequences; i++)
    {
      seqinfo_s * p = d->seqindex + i;
      p->seq = r;
      r += p->seqlen;
      progress_update(i+1);
    }
  progress_done();
}

void db_hash(struct db * d)
{
  progress_init("Computing hashes: ", d->sequences);
  for(uint64_t i = 0; i < d->sequences; i++)
    {
      seqinfo_s * p = d->seqindex + i;
      d->seqindex[i].hash = zobrist_hash((unsigned char *)(p->seq),
                                         p->seqlen,
                                         p->v_gene_no,
                                         p->j_gene_no);
      progress_update(i+1);
    }
  progress_done();
}

void db_free(struct db * d)
{
  if (d->residues_p)
    xfree(d->residues_p);
  if (d->seqindex)
    {
      for (uint64_t i = 0; i < d->sequences; i++)
        {
          if (d->seqindex[i].sequence_id)
            {
              xfree(d->seqindex[i].sequence_id);
              d->seqindex[i].sequence_id = nullptr;
            }
          if (d->seqindex[i].keep)
            {
              xfree(d->seqindex[i].keep);
              d->seqindex[i].keep = nullptr;
            }
        }
      xfree(d->seqindex);
    }
  d->repertoire_id_vector.clear();
  d->repertoire_id_map.clear();
  delete d;
}

uint64_t db_getsequencecount(struct db * d)
{
  return d->sequences;
}

uint64_t db_getresiduescount(struct db * d)
{
  return d->residues_count;
}

unsigned int db_getlongestsequence(struct db * d)
{
  return d->longest;
}

uint64_t db_get_repertoire_count(struct db * d)
{
  return d->repertoire_count;
}

uint64_t db_gethash(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].hash;
}

char * db_getsequence(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].seq;
}

unsigned int db_getsequencelen(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].seqlen;
}

uint64_t db_get_v_gene(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].v_gene_no;
}

uint64_t db_get_j_gene(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].j_gene_no;
}

uint64_t db_get_count(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].count;
}

int db_get_repertoire_id_no(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].repertoire_id_no;
}

const char * db_get_repertoire_id(struct db * d, int repertoire_id_no)
{
  return d->repertoire_id_vector[repertoire_id_no].c_str();
}

char * db_get_sequence_id(struct db * d, uint64_t seqno)
{
  char * sid = d->seqindex[seqno].sequence_id;
  if (sid)
    return sid;
  else
    return (char *) EMPTYSTRING;
}

uint64_t db_get_v_gene_count()
{
  return v_gene_vector.size();
}

uint64_t db_get_j_gene_count()
{
  return j_gene_vector.size();
}

const char * db_get_v_gene_name(struct db * d, uint64_t seqno)
{
  int v_gene_no = d->seqindex[seqno].v_gene_no;
  return v_gene_vector[v_gene_no].c_str();
}

const char * db_get_j_gene_name(struct db * d, uint64_t seqno)
{
  int j_gene_no = d->seqindex[seqno].j_gene_no;
  return j_gene_vector[j_gene_no].c_str();
}

void db_fprint_sequence(FILE * f, struct db * d, uint64_t seqno)
{
  char * seq = db_getsequence(d, seqno);
  unsigned int len = db_getsequencelen(d, seqno);
  if (opt_nucleotides)
    {
      for (unsigned int i = 0; i < len; i++)
        fputc(nt_chars[(int)(seq[i])], f);
    }
  else
    {
      for (unsigned int i = 0; i < len; i++)
        fputc(aa_chars[(int)(seq[i])], f);
    }
}

char * db_get_keep_columns(struct db * d, uint64_t seqno)
{
  char * keep = d->seqindex[seqno].keep;
  if (keep)
    return keep;
  else
    return (char *) EMPTYSTRING;
}
