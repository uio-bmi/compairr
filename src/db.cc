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
  std::vector<std::string> repertoire_id_vector;
  std::map<std::string, int> repertoire_id_map;
  int col_junction;
  int col_junction_aa;
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
  d->col_duplicate_count = 0;
  d->col_v_call = 0;
  d->col_j_call = 0;
  d->col_repertoire_id = 0;
  d->col_sequence_id = 0;

  return d;
}

void parse_airr_tsv_header(char * line,
                           struct db * d,
                           bool require_repertoire_id,
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
      i++;
    }

  if (! (d->col_repertoire_id   || ! require_repertoire_id) ||
      ! (d->col_sequence_id     || ! require_sequence_id)   ||
      ! (d->col_duplicate_count ||   opt_ignore_counts)     ||
      ! (d->col_v_call          ||   opt_ignore_genes)      ||
      ! (d->col_j_call          ||   opt_ignore_genes)      ||
      ! (d->col_junction        || ! opt_nucleotides)       ||
      ! (d->col_junction_aa     ||   opt_nucleotides))
    {
      fprintf(logfile,
        "\nMissing essential column(s) in header of AIRR TSV input file:");

      if (! (d->col_repertoire_id   || ! require_repertoire_id))
        fprintf(logfile, " repertoire_id");
      if (! (d->col_sequence_id     || ! require_sequence_id))
        fprintf(logfile, " sequence_id");
      if (! (d->col_duplicate_count || opt_ignore_counts))
        fprintf(logfile, " duplicate_count");
      if (! (d->col_v_call          || opt_ignore_genes))
        fprintf(logfile, " v_call");
      if (! (d->col_j_call          || opt_ignore_genes))
        fprintf(logfile, " j_call");
      if (! (d->col_junction        || ! opt_nucleotides))
        fprintf(logfile, " junction");
      if (! (d->col_junction_aa     || opt_nucleotides))
        fprintf(logfile, " junction_aa");

      fprintf(logfile, "\n");
      exit(1);
    }
}

void parse_airr_tsv_line(char * line,
                         uint64_t lineno,
                         struct db * d,
                         bool require_repertoire_id,
                         bool require_sequence_id)
{
  char * repertoire_id = nullptr;
  char * sequence_id = nullptr;
  char * duplicate_count = nullptr;
  char * v_call = nullptr;
  char * j_call = nullptr;
  char * junction = nullptr;
  char * junction_aa = nullptr;

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


  /* handle repertoire_id */

  if (repertoire_id && *repertoire_id)
    {
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
    }
  else if (require_repertoire_id)
    {
      fprintf(logfile,
              "\n\nError: missing or empty repertoire_id value on line: %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }
  else
    {
      p->repertoire_id_no = 0;
    }


  /* handle sequence_id */

  if (sequence_id && *sequence_id)
    {
      p->sequence_id = xstrdup(sequence_id);
    }
  else if (require_sequence_id)
    {
      fprintf(logfile,
              "\n\nError: missing or empty sequence_id value on line: %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }
  else
    {
      p->sequence_id = nullptr;
    }


  /* handle duplicate_count */

  long count = 1;

  if (duplicate_count && *duplicate_count)
    {
      char * endptr;
      count = strtol(duplicate_count, &endptr, 10);
      if ((endptr == nullptr) || (*endptr) || (count < 1))
        {
          fprintf(logfile, "\n\nError: Illegal duplicate_count on line %"
                  PRIu64 ": %s\n", lineno, duplicate_count);
          exit(1);
        }
    }
  else if (! opt_ignore_counts)
    {
      fprintf(logfile,
              "\n\nError: missing or empty duplicate_count on line: %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }

  p->count = count;
  d->total_duplicate_count += count;


  /* handle v_call */

  if (v_call && *v_call)
    {
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
    }
  else if (opt_ignore_genes)
    {
      p->v_gene_no = 0;
    }
  else
    {
      fprintf(logfile,
              "\n\nError: missing or empty v_call value on line: %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }


  /* handle j_call */

  if (j_call && *j_call)
    {
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
    }
  else if (opt_ignore_genes)
    {
      p->j_gene_no = 0;
    }
  else
    {
      fprintf(logfile,
              "\n\nError: missing or empty j_call value on line: %"
              PRIu64 "\n",
              lineno);
      exit(1);
    }


  /* handle junction and junction_aa */

  if (opt_nucleotides)
    {
      if (! (junction && *junction))
        {
          fprintf(logfile,
                  "\n\nError: missing or empty junction value on line: %"
                  PRIu64 "\n",
                  lineno);
          exit(1);
        }

    }
  else
    {
      if (! (junction_aa && *junction_aa))
        {
          fprintf(logfile,
                  "\n\nError: missing or empty junction_aa value on line: %"
                  PRIu64 "\n",
                  lineno);
          exit(1);
        }
    }

  /* make room for more residues */

  unsigned int len_estimate = 0;
  if (opt_nucleotides)
    len_estimate = strlen(junction);
  else
    len_estimate = strlen(junction_aa);

  if (d->residues_count + len_estimate > d->residues_alloc)
    {
      d->residues_alloc += MEMCHUNK;
      d->residues_p = static_cast<char *>
        (xrealloc(d->residues_p, d->residues_alloc));
    }

  /* scan and store sequence */

  char * q = d->residues_p + d->residues_count;
  unsigned int seqlen = 0;

  if (opt_nucleotides)
    {
      for(unsigned int i = 0; i < len_estimate; i++)
        {
          unsigned char c = junction[i];
          signed char m = map_nt[static_cast<unsigned int>(c)];
          if (m >= 0)
            {
              *q++ = m;
              seqlen++;
            }
          else
            {
              if ((c >= 32) && (c <= 126))
                fprintf(logfile,
                        "\n\nError: Illegal character '%c' in sequence "
                        "on line %" PRIu64 "\n",
                        c,
                        lineno);
              else
                fprintf(logfile,
                        "\n\nError: Illegal character (ascii no %d) in sequence "
                        "on line %" PRIu64 "\n",
                        c,
                        lineno);
              exit(1);
            }
        }
    }
  else
    {
      for(unsigned int i = 0; i < len_estimate; i++)
        {
          unsigned char c = junction_aa[i];
          signed char m = map_aa[static_cast<unsigned int>(c)];
          if (m >= 0)
            {
              *q++ = m;
              seqlen++;
            }
          else
            {
              if ((c >= 32) && (c <= 126))
                fprintf(logfile,
                        "\n\nError: Illegal character '%c' in sequence "
                        "on line %" PRIu64 "\n",
                        c,
                        lineno);
              else
                fprintf(logfile,
                        "\n\nError: Illegal character (ascii no %d) in sequence "
                        "on line %" PRIu64 "\n",
                        c,
                        lineno);
              exit(1);
            }
        }
    }

  d->residues_count += seqlen;
  p->seqlen = seqlen;
  if (seqlen > d->longest)
    d->longest = seqlen;
  if (seqlen < d->shortest)
    d->shortest = seqlen;

  p->hash = 0;

  d->sequences++;
}

void db_read(struct db * d,
             const char * filename,
             bool require_repertoire_id,
             bool require_sequence_id)
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
                                    require_repertoire_id,
                                    require_sequence_id);
              state = 1;
            }
        }
      else
        {
          parse_airr_tsv_line(line,
                              lineno,
                              d,
                              require_repertoire_id,
                              require_sequence_id);
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

  /* if repertoire_id is not used, insert empty string as ID of first */

  if (d->repertoire_count == 0)
    {
      d->repertoire_id_vector.push_back("");
      d->repertoire_id_map.insert({"", 0});
      d->repertoire_count = 1;
    }

  fprintf(logfile,
          "Repertoires:       %" PRIu64 "\n"
          "Sequences:         %" PRIu64 "\n"
          "Residues:          %" PRIu64 "\n"
          "Shortest:          %u\n"
          "Longest:           %u\n"
          "Average length:    %4.1lf\n"
          "Total dupl. count: %" PRIu64 "\n",
          d->repertoire_count,
          d->sequences,
          d->residues_count,
          d->shortest,
          d->longest,
          1.0 * d->residues_count / d->sequences,
          d->total_duplicate_count);

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
    xfree(d->seqindex);
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
  return v_gene_vector[d->seqindex[seqno].v_gene_no].c_str();
}

const char * db_get_j_gene_name(struct db * d, uint64_t seqno)
{
  return j_gene_vector[d->seqindex[seqno].j_gene_no].c_str();
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
