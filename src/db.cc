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

#include "vdjsearch.h"

/* How much memory for residues and sequences should we allocate each time? */

#define MEMCHUNK 1048576
#define SEQCHUNK 65536

constexpr int LINEALLOC = 2048;

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

struct seqinfo_s
{
  uint64_t hash;
  char * seq;
  short seqlen;
  unsigned int sample_no;
  short v_gene_no;
  short d_gene_no;
  double freq;
};

typedef struct seqinfo_s seqinfo_t;

struct db
{
  seqinfo_t * seqindex;
  uint64_t seqindex_alloc;
  unsigned int sequences;
  unsigned int longest;
  unsigned int shortest;
  char * residues_p;
  uint64_t residues_alloc;
  uint64_t residues_count;
  char * * sample_list;
  uint64_t sample_alloc;
  uint64_t sample_count;
};

static char * * v_gene_list = nullptr;
static uint64_t v_gene_alloc = 0;
static uint64_t v_gene_count = 0;
static char * * d_gene_list = nullptr;
static uint64_t d_gene_alloc = 0;
static uint64_t d_gene_count = 0;

uint64_t list_insert(char * * * list, uint64_t * alloc, uint64_t * count, char * item)
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

int compare_bysample(const void * a, const void * b)
{
  const seqinfo_s * x = (const seqinfo_s *) a;
  const seqinfo_s * y = (const seqinfo_s *) b;

  if (x->sample_no < y->sample_no)
    return -1;
  else if (x->sample_no > y->sample_no)
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
  d->sample_list = nullptr;
  d->sample_alloc = 0;
  d->sample_count = 0;

  return d;
}

void db_read(struct db * d, const char * filename)
{
  FILE * fp = nullptr;
  if (filename)
    {
      fp = fopen_input(filename);
      if (!fp)
        {
          fprintf(stderr, "\nError: Unable to open input data file (%s).\n", filename);
          exit(1);
        }
    }
  else
    fp = stdin;

  /* get file size */

  struct stat fs;

  if (fstat(fileno(fp), & fs))
    {
      fprintf(stderr, "\nUnable to fstat on input file (%s)\n", filename);
      exit(1);
    }
  bool is_regular = S_ISREG(fs.st_mode);
  int64_t filesize = is_regular ? fs.st_size : 0;

  if (! is_regular)
    fprintf(logfile, "Waiting for data... (Hit Ctrl-C and run %s -h if you meant to read data from a file.)\n", PROG_NAME);

  char line[LINEALLOC];
  line[0] = 0;
  if (!fgets(line, LINEALLOC, fp))
    line[0] = 0;

  unsigned int lineno = 1;

  progress_init("Reading sequences:", static_cast<uint64_t>(filesize));

  d->longest = 0;
  d->shortest = UINT_MAX;

  while(line[0])
    {
      /* make room for another sequence */

      if (d->sequences >= d->seqindex_alloc)
        {
          d->seqindex_alloc += SEQCHUNK;
          d->seqindex = static_cast<seqinfo_t *>
            (xrealloc(d->seqindex, d->seqindex_alloc * sizeof(seqinfo_s)));
        }

      char * aa_seq = NULL;
      char * freq_s = NULL;
      char * v_gene = NULL;
      char * d_gene = NULL;
      char * sample = NULL;

      double freq = 0.0;

      char tab[] = "\t";
      char nl[] = "\r\n";

      aa_seq = strtok(line, tab);

      if (aa_seq)
        freq_s = strtok(NULL, tab);

      if (freq_s)
        {
          freq = strtod(freq_s, NULL);
          v_gene = strtok(NULL, tab);
        }

      if (v_gene)
        d_gene = strtok(NULL, tab);

      if (d_gene)
        sample = strtok(NULL, nl);

      if (! (aa_seq && freq_s && v_gene && d_gene && sample))
        {
          fprintf(stderr, "Missing data on line: %d\n", lineno);
          fprintf(stderr, "aa_seq: %s freq: %s v_gene: %s d_gene: %s sample: %s\n",
                 aa_seq, freq_s, v_gene, d_gene, sample);
          fatal("fatal");
        }

      /* make room for more residues */

      unsigned int len_estimate = freq_s - aa_seq - 1;

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
          unsigned char c = aa_seq[i];
          signed char m = map_aa[static_cast<unsigned int>(c)];
          if (m >= 0)
            {
              *q++ = m;
              seqlen++;
            }
          else
            {
              if ((c >= 32) && (c <= 126))
                fprintf(stderr,
                        "\nError: Illegal character '%c' in sequence on line %u\n",
                        c,
                        lineno);
              else
                fprintf(stderr,
                        "\nError: Illegal character (ascii no %d) in sequence on line %u\n",
                        c,
                        lineno);
              exit(1);
            }
        }

      uint64_t v_gene_index = list_insert(& v_gene_list,
                                          & v_gene_alloc,
                                          & v_gene_count,
                                          v_gene);

      uint64_t d_gene_index = list_insert(& d_gene_list,
                                          & d_gene_alloc,
                                          & d_gene_count,
                                          d_gene);

      uint64_t sample_index = list_insert(& d->sample_list,
                                          & d->sample_alloc,
                                          & d->sample_count,
                                          sample);

      p->seqlen = seqlen;
      p->sample_no = sample_index;
      p->v_gene_no = v_gene_index;
      p->d_gene_no = d_gene_index;
      p->freq = freq;
      p->hash = 0;

      if (seqlen > d->longest)
        d->longest = seqlen;
      if (seqlen < d->shortest)
        d->shortest = seqlen;

      d->sequences++;
      d->residues_count += seqlen;

      /* get next line */

      line[0] = 0;
      if (!fgets(line, LINEALLOC, fp))
        line[0] = 0;
      lineno++;

      if (is_regular)
        progress_update(static_cast<uint64_t>(ftell(fp)));
    }
  progress_done();

  fclose(fp);

  fprintf(logfile,
          "Sequences: %u, residues: %llu, shortest: %u, longest: %u, average: %4.1lf\n",
          d->sequences,
          d->residues_count,
          d->shortest,
          d->longest,
          1.0 * d->residues_count / d->sequences);

  fprintf(logfile, "Samples:           %llu\n", d->sample_count);

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

  /* sort sequences by sample */

  progress_init("Sorting:          ", 1);
  qsort(d->seqindex, d->sequences, sizeof(seqinfo_s), compare_bysample);
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
                                         p->d_gene_no);
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
  xfree(d);
}

unsigned int db_getsequencecount(struct db * d)
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

uint64_t db_getsamplecount(struct db * d)
{
  return d->sample_count;
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

uint64_t db_get_d_gene(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].d_gene_no;
}

double db_get_freq(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].freq;
}

uint64_t db_getsampleno(struct db * d, uint64_t seqno)
{
  return d->seqindex[seqno].sample_no;
}

char * db_getsamplename(struct db * d, uint64_t sample_no)
{
  return d->sample_list[sample_no];
}

uint64_t db_get_v_gene_count()
{
  return v_gene_count;
}

uint64_t db_get_d_gene_count()
{
  return d_gene_count;
}
