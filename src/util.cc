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

static const char * progress_prompt;
static uint64_t progress_next;
static uint64_t progress_size;
static uint64_t progress_chunk;
static const uint64_t progress_granularity = 200;
const size_t memalignment = 16;
static std::chrono::time_point<std::chrono::steady_clock> time_point_start;

void progress_init(const char * prompt, uint64_t size)
{
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ?
    1 : size / progress_granularity;
  progress_next = progress_chunk;
  if (opt_log)
    fprintf(logfile, "%s", prompt);
  else
    fprintf(logfile, "%s %.0f%%", prompt, 0.0);
  fflush(logfile);
  time_point_start = std::chrono::steady_clock::now();
}

void progress_update(uint64_t progress)
{
  if ((!opt_log) && (progress >= progress_next))
    {
      fprintf(logfile, "  \r%s %.0f%%", progress_prompt,
              100.0 * static_cast<double>(progress)
              / static_cast<double>(progress_size));
      progress_next = progress + progress_chunk;
      fflush(logfile);
    }
}

void progress_done()
{
  auto time_point_now = std::chrono::steady_clock::now();
  double time_diff = 0.000000001 * (time_point_now - time_point_start)
    / std::chrono::nanoseconds(1);
  if (opt_log)
    fprintf(logfile, " %.0f%% (%.9lfs)\n", 100.0, time_diff);
  else
    fprintf(logfile, "  \r%s %.0f%% (%.9lfs)\n", progress_prompt, 100.0,
            time_diff);
  fflush(logfile);
}

int64_t gcd(int64_t a, int64_t b)
{
  if (b == 0)
  {
    return a;
  }
  else
  {
    return gcd(b, a % b);
  }
}

[[ noreturn ]] void fatal(const char * msg)
{
  fprintf(stderr, "\nError: %s\n", msg);
  exit(1);
}

void * xmalloc(size_t size)
{
  if (size == 0)
    size = 1;
  void * t = nullptr;
#ifdef _WIN32
  t = _aligned_malloc(size, memalignment);
#else
  if (posix_memalign(& t, memalignment, size))
    t = nullptr;
#endif
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  if (size == 0)
    size = 1;
#ifdef _WIN32
  void * t = _aligned_realloc(ptr, size, memalignment);
#else
  void * t = realloc(ptr, size);
#endif
  if (!t)
    fatal("Unable to reallocate enough memory.");
  return t;
}

void xfree(void * ptr)
{
  if (ptr)
    {
#ifdef _WIN32
      _aligned_free(ptr);
#else
      free(ptr);
#endif
    }
  else
    fatal("Trying to free a null pointer");
}

char * xstrdup(const char * s)
{
  char * t = strdup(s);
  if (t == nullptr)
    fatal("Out of memory");
  return t;
}

FILE * fopen_input(const char * filename)
{
  /* open the input stream given by filename, but use stdin if name is - */
  if (strcmp(filename, "-") == 0)
    {
      int fd = dup(STDIN_FILENO);
      if (fd < 0)
        return nullptr;
      else
        return fdopen(fd, "rb");
    }
  else
    return fopen(filename, "rb");
}

FILE * fopen_output(const char * filename)
{
  /* open the output stream given by filename, but use stdout if name is - */
  if (strcmp(filename, "-") == 0)
    {
      int fd = dup(STDOUT_FILENO);
      if (fd < 0)
        return nullptr;
      else
        return fdopen(fd, "w");
    }
  else
    return fopen(filename, "w");
}

int64_t seq_diff(unsigned char * a, unsigned char * b, int64_t len)
{
  /* Count number of different characters in a and b of length len */
  int64_t diffs = 0;
  for (int64_t i = 0; i < len; i++)
    if (*a++ != *b++)
      {
        diffs++;
        if (diffs > opt_differences)
          break;
      }
  return diffs;
}
