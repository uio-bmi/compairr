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

uint64_t arch_get_memused()
{
#ifdef _WIN32

  PROCESS_MEMORY_COUNTERS pmc;
  GetProcessMemoryInfo(GetCurrentProcess(),
                       &pmc,
                       sizeof(PROCESS_MEMORY_COUNTERS));
  return pmc.PeakWorkingSetSize;

#else

  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);

# ifdef __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return static_cast<uint64_t>(r_usage.ru_maxrss);
# else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return static_cast<uint64_t>(r_usage.ru_maxrss * 1024);
# endif

#endif
}

uint64_t arch_get_memtotal()
{
#ifdef _WIN32

  MEMORYSTATUSEX ms;
  ms.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&ms);
  return ms.ullTotalPhys;

#elif defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  size_t length = sizeof(ram);
  if(sysctl(mib, 2, &ram, &length, nullptr, 0) == -1)
    fatal("Cannot determine amount of RAM");
  return static_cast<uint64_t>(ram);

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  int64_t phys_pages = sysconf(_SC_PHYS_PAGES);
  int64_t pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) || (pagesize == -1))
    fatal("Cannot determine amount of RAM");
  return static_cast<uint64_t>(pagesize * phys_pages);

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}

void arch_srandom(unsigned int seed)
{
  /* initialize pseudo-random number generator */

#ifdef _WIN32
  srand(seed);
#else
  srandom(seed);
#endif
}

uint64_t arch_random()
{
#ifdef _WIN32
  return static_cast<uint64_t>(rand());
#else
  return static_cast<uint64_t>(random());
#endif
}
