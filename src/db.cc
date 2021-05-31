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

const static char * aa_chars = "ACDEFGHIKLMNPQRSTVWY";

struct seqinfo_s
{
  uint64_t hash;
  char * seq;
  short seqlen;
  unsigned int repertoire_id_no;
  short v_gene_no;
  short j_gene_no;
  uint64_t count;
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
  char * * repertoire_id_list;
  uint64_t repertoire_id_alloc;
  uint64_t repertoire_count;
};

static char * * v_gene_list = 0;
static uint64_t v_gene_alloc = 0;
static uint64_t v_gene_count = 0;
static char * * j_gene_list = 0;
static uint64_t j_gene_alloc = 0;
static uint64_t j_gene_count = 0;

void db_init()
{
  v_gene_list = 0;
  v_gene_alloc = 0;
  v_gene_count = 0;
  j_gene_list = 0;
  j_gene_alloc = 0;
  j_gene_count = 0;
}

void db_exit()
{
  if (v_gene_list)
    {
      for (uint64_t i = 0; i < v_gene_count; i++)
	{
	  xfree(v_gene_list[i]);
	  v_gene_list[i] = 0;
	}
      xfree(v_gene_list);
      v_gene_list = 0;
    }

  if (j_gene_list)
    {
      for (uint64_t i = 0; i < j_gene_count; i++)
	{
	  xfree(j_gene_list[i]);
	  j_gene_list[i] = 0;
	}
      xfree(j_gene_list);
      j_gene_list = 0;
    }
}

uint64_t list_insert(char * * * list,
                     uint64_t * alloc,
                     uint64_t * count,
                     char * item)
{
  /* linear list */

  /* try to find it */
  uint64_t c = *count;
  for(uint64_t i = 0; i < c; i++)
    if (strcmp(item, (*list)[i]) == 0)
        return i;

  /* expand if necessary */
  if (c >= *alloc)
    {
      *alloc += 256;
      *list = static_cast<char **>(xrealloc(*list, *alloc * sizeof(char*)));
    }

  /* insert */
  (*list)[c] = strdup(item);
  (*count)++;
  return c;
}

int compare_byrepertoire_id(const void * a, const void * b)
{
  const seqinfo_s * x = (const seqinfo_s *) a;
  const seqinfo_s * y = (const seqinfo_s *) b;

  if (x->repertoire_id_no < y->repertoire_id_no)
    return -1;
  else if (x->repertoire_id_no > y->repertoire_id_no)
    return +1;
  else
    return 0;
}

struct db * db_create()
{
  struct db * d = (struct db *) xmalloc(sizeof(struct db));

  d->seqindex = nullptr;
  d->seqindex_alloc = 0;
  d->sequences = 0;
  d->longest = 0;
  d->shortest = UINT_MAX;
  d->residues_p = nullptr;
  d->residues_alloc = 0;
  d->residues_count = 0;
  d->repertoire_id_list = nullptr;
  d->repertoire_id_alloc = 0;
  d->repertoire_count = 0;

  return d;
}

static int col_junction_aa = 0;
static int col_duplicate_count = 0;
static int col_v_call = 0;
static int col_j_call = 0;
static int col_repertoire_id = 0;

void parse_airr_tsv_header(char * line)
{
  char delim[] = "\t";
  char * string = line;
  char * token = nullptr;

  int i = 1;

  while ((token = strsep(& string, delim)) != NULL)
    {
      if (strcmp(token, "junction_aa") == 0)
        {
          col_junction_aa = i;
        }
      else if (strcmp(token, "duplicate_count") == 0)
        {
          col_duplicate_count = i;
        }
      else if (strcmp(token, "v_call") == 0)
        {
          col_v_call = i;
        }
      else if (strcmp(token, "j_call") == 0)
        {
          col_j_call = i;
        }
      else if (strcmp(token, "repertoire_id") == 0)
        {
          col_repertoire_id = i;
        }
      i++;
    }

  if (! col_junction_aa ||
      ! col_duplicate_count ||
      ! col_v_call ||
      ! col_j_call ||
      ! col_repertoire_id)
    {
      fprintf(logfile,
        "\nMissing essential column(s) in header of AIRR TSV input file:");
      if (! col_junction_aa)
        fprintf(logfile, " junction_aa");
      if (! col_duplicate_count)
        fprintf(logfile, " duplicate_count");
      if (! col_v_call)
        fprintf(logfile, " v_call");
      if (! col_j_call)
        fprintf(logfile, " j_call");
      if (! col_repertoire_id)
        fprintf(logfile, " repertoire_id");
      fprintf(logfile, "\n");
      exit(1);
    }
}

void parse_airr_tsv_line(char * line, uint64_t lineno, struct db * d)
{
  char * junction_aa = NULL;
  char * duplicate_count = NULL;
  char * v_call = NULL;
  char * j_call = NULL;
  char * repertoire_id = NULL;

  char delim[] = "\t";
  char * string = line;
  char * token = nullptr;

  int i = 1;

  while ((token = strsep(& string, delim)) != NULL)
    {
      if (i == col_junction_aa)
        {
          junction_aa = token;
        }
      else if (i == col_duplicate_count)
        {
          duplicate_count = token;
        }
      else if (i == col_v_call)
        {
          v_call = token;
        }
      else if (i == col_j_call)
        {
          j_call = token;
        }
      else if (i == col_repertoire_id)
        {
          repertoire_id = token;
        }
      i++;
    }

  /* check that all values are read */

  if (! (junction_aa && duplicate_count && v_call && j_call && repertoire_id))
    {
      fprintf(logfile,
              "\n\nError: Missing data on line: %" PRIu64 "\n",
              lineno);
      fprintf(logfile,
              "junction_aa: %s duplicate_count: %s "
              "v_call: %s j_call: %s repertoire_id: %s\n",
              junction_aa, duplicate_count, v_call, j_call, repertoire_id);
      exit(1);
    }

  /* check that strings are not empty */

  if (! (*junction_aa && *duplicate_count && *v_call && *j_call &&
        *repertoire_id))
    {
      fprintf(logfile, "\n\nError: Empty string on line: %" PRIu64 "\n",
              lineno);
      fprintf(logfile, "junction_aa: %s duplicate_count: %s "
              "v_call: %s j_call: %s repertoire_id: %s\n",
              junction_aa, duplicate_count, v_call, j_call, repertoire_id);
      exit(1);
    }

  /* check the converted duplicate_count */

  long count = 0;
  char * endptr;
  count = strtol(duplicate_count, &endptr, 10);

  if ((endptr == duplicate_count) || (endptr == nullptr) || (*endptr) ||
      (count <= 0))
    {
      fprintf(logfile, "\n\nError: Illegal duplicate_count on line %"
              PRIu64 ": %s\n", lineno, duplicate_count);
      exit(1);
    }

  /* make room for another sequence */

  if (d->sequences >= d->seqindex_alloc)
    {
      d->seqindex_alloc += SEQCHUNK;
      d->seqindex = static_cast<seqinfo_t *>
        (xrealloc(d->seqindex, d->seqindex_alloc * sizeof(seqinfo_s)));
    }

  /* make room for more residues */

  unsigned int len_estimate = strlen(junction_aa);

  if (d->residues_count + len_estimate > d->residues_alloc)
    {
      d->residues_alloc += MEMCHUNK;
      d->residues_p = static_cast<char *>
        (xrealloc(d->residues_p, d->residues_alloc));
    }

  /* scan and store sequence */

  seqinfo_t * p = d->seqindex + d->sequences;
  char * q = d->residues_p + d->residues_count;
  unsigned int seqlen = 0;

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

  uint64_t v_gene_index = list_insert(& v_gene_list,
                                      & v_gene_alloc,
                                      & v_gene_count,
                                      v_call);

  uint64_t j_gene_index = list_insert(& j_gene_list,
                                      & j_gene_alloc,
                                      & j_gene_count,
                                      j_call);

  uint64_t repertoire_id_index = list_insert(& d->repertoire_id_list,
                                      & d->repertoire_id_alloc,
                                      & d->repertoire_count,
                                      repertoire_id);

  p->seqlen = seqlen;
  p->repertoire_id_no = repertoire_id_index;
  p->v_gene_no = v_gene_index;
  p->j_gene_no = j_gene_index;
  p->count = count;
  p->hash = 0;

  if (seqlen > d->longest)
    d->longest = seqlen;
  if (seqlen < d->shortest)
    d->shortest = seqlen;

  d->sequences++;
  d->residues_count += seqlen;
}

void db_read(struct db * d, const char * filename)
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
  int64_t filesize = is_regular ? fs.st_size : 0;

  if (! is_regular)
    fprintf(logfile, "Waiting for data from standard input...\n");

  size_t line_alloc = 4096;
  char * line = (char *) xmalloc(line_alloc);
  uint64_t lineno = 0;
  ssize_t linelen = 0;

  d->longest = 0;
  d->shortest = UINT_MAX;

  int state = 0;

  progress_init("Reading sequences:", static_cast<uint64_t>(filesize));

  linelen = getline(& line, & line_alloc, fp);

  if (linelen < 0)
    fatal("Unable to read from the input file");

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
              parse_airr_tsv_header(line);
              state = 1;
            }
        }
      else
        {
          parse_airr_tsv_line(line, lineno, d);
        }

      /* update progress */

      if (is_regular)
        progress_update(static_cast<uint64_t>(ftell(fp)));

      /* get next line */

      linelen = getline(& line, & line_alloc, fp);

      if (linelen < 0)
        break;

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

  fprintf(logfile,
          "Sequences:         %" PRIu64 "\n"
          "Residues:          %" PRIu64 "\n"
          "Shortest:          %u\n"
          "Longest:           %u\n"
          "Average length:    %4.1lf\n",
          d->sequences,
          d->residues_count,
          d->shortest,
          d->longest,
          1.0 * d->residues_count / d->sequences);

  fprintf(logfile, "Repertoires:       %" PRIu64 "\n", d->repertoire_count);

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

  /* sort sequences by repertoire_id */

  progress_init("Sorting:          ", 1);
  qsort(d->seqindex, d->sequences, sizeof(seqinfo_s), compare_byrepertoire_id);
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
  if (d->repertoire_id_list)
    {
      for (uint64_t i = 0; i < d->repertoire_count; i++)
	{
	  xfree(d->repertoire_id_list[i]);
	  d->repertoire_id_list[i] = 0;
	}
      xfree(d->repertoire_id_list);
      d->repertoire_id_list = 0;
    }

  if (d->residues_p)
    xfree(d->residues_p);
  if (d->seqindex)
    xfree(d->seqindex);
  xfree(d);
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
  return opt_ignore_counts ? 1 : d->seqindex[seqno].count;
}

uint64_t db_get_repertoire_id_no(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].repertoire_id_no;
}

char * db_get_repertoire_id(struct db * d, uint64_t repertoire_id_no)
{
  return d->repertoire_id_list[repertoire_id_no];
}

uint64_t db_get_v_gene_count()
{
  return v_gene_count;
}

uint64_t db_get_j_gene_count()
{
  return j_gene_count;
}

char * db_get_v_gene_name(struct db * d, uint64_t seqno)
{
  return v_gene_list[d->seqindex[seqno].v_gene_no];
}

char * db_get_j_gene_name(struct db * d, uint64_t seqno)
{
  return j_gene_list[d->seqindex[seqno].j_gene_no];
}

void db_fprint_sequence(FILE * f, struct db * d, uint64_t seqno)
{
  char * seq = db_getsequence(d, seqno);
  unsigned int len = db_getsequencelen(d, seqno);
  for (unsigned int i = 0; i < len; i++)
    fputc(aa_chars[(int)(seq[i])], f);
}
