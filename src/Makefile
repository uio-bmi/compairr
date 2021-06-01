#    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
#    Department of Informatics, University of Oslo,
#    PO Box 1080 Blindern, NO-0316 Oslo, Norway

# Makefile for CompAIRR

COMMON = -g -std=c++11
EXTRAOBJ =
LINKOPT =
LIBS = -lpthread
WARNINGS = -Wall -Wextra \
#-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic

# Run "make RELEASE=1" to compile for release
ifdef RELEASE
	COMMON += -DNDEBUG
endif

# Run "make PROFILE=1" to compile for profiling
ifdef PROFILE
	COMMON += -pg
endif

# Run "make COVERAGE=1" to compile for coverage tests
ifdef COVERAGE
	COMMON += -fprofile-arcs -ftest-coverage -O0
	LIBS += -lgcov
	LINKFLAGS += --coverage
else
	COMMON += -flto -O3
endif

# Machine specific
MACHINE = $(shell uname -m)
ifeq ($(MACHINE), x86_64)
	COMMON += -march=x86-64 -mtune=generic
else ifeq ($(MACHINE), arm64)
	COMMON += -march=armv8-a+simd -mtune=generic
else ifeq ($(MACHINE), aarch64)
	COMMON += -march=armv8-a+simd -mtune=generic
else ifeq ($(MACHINE), ppc64le)
	COMMON += -mcpu=power8
endif

# OS specific
ifeq ($(CXX), x86_64-w64-mingw32-g++)
	LIBS += -lpsapi
	LINKOPT += -static
else
	WARNINGS += -pedantic
endif

LINKFLAGS = $(COMMON) $(LINKOPT)

CXXFLAGS = $(COMMON) $(WARNINGS)

PROG = compairr

OBJS = arch.o bloompat.o cluster.o compairr.o db.o hashtable.o \
	overlap.o util.o variants.o zobrist.o

DEPS = Makefile threads.h \
	arch.h bloompat.h cluster.h compairr.h db.h hashtable.h \
	overlap.h util.h variants.h zobrist.h

all : $(PROG)

compairr : $(OBJS) $(DEPS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean :
	rm -rf compairr *.o *~ gmon.out *.gcno *.gcda *.gcov

.o : .cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
