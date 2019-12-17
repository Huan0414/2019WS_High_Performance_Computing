# Basic Defintions for using GNU-compiler suite sequentially
# requires setting of COMPILER=CLANG_

CC	= clang
CXX     = clang++
#CXX     = /usr/bin/clang++-6.0
#F77	= gfortran
LINKER  = ${CXX}



#http://clang.llvm.org/docs/UsersManual.html#options-to-control-error-and-warning-messages
#WARNINGS = -Weverything -Wno-c++98-compat -ferror-limit=3
WARNINGS =  -Weverything -Wno-sign-conversion -Wno-c++98-compat -Wno-c++98-compat-bind-to-temporary-copy -Wno-shorten-64-to-32 -Wno-padded
#-fsyntax-only -Wdocumentation -Wconversion -Wshadow -Wfloat-conversion -pedantic
CXXFLAGS += -Ofast -std=c++17 -ferror-limit=1 ${WARNINGS}
# -ftrapv

LINKFLAGS +=

#   very good check
# http://clang.llvm.org/extra/clang-tidy/
#TIDYFLAGS = -checks='modernize*'
#   good check, see:  http://llvm.org/docs/CodingStandards.html#include-style
TIDYFLAGS = -checks=llvm*,-llvm-header-guard -header-filter=.* -enable-check-profile -extra-arg="-std=c++17"
#TIDYFLAGS = -checks=llvm*,readability-*,-llvm-header-guard  -header-filter=.* -export-fixes=fixes.txt
#
#TIDYFLAGS = -checks='readability-*'  -header-filter=.*
#   ???
#TIDYFLAGS = -checks='cert*'  -header-filter=.*
#   MPI checks ??
#TIDYFLAGS = -checks='mpi*'
#   ??
#TIDYFLAGS = -checks='performance*'   -header-filter=.*
#TIDYFLAGS = -checks='portability-*'  -header-filter=.*
#TIDYFLAGS = -checks='readability-*'  -header-filter=.*

default: ${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@

clean:
	@rm -f ${PROGRAM} ${OBJECTS}

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar

tidy_check:
	clang-tidy ${SOURCES} ${TIDYFLAGS} -- ${SOURCES}
# see also http://clang-developers.42468.n3.nabble.com/Error-while-trying-to-load-a-compilation-database-td4049722.html


run: clean ${PROGRAM}
#	time  ./${PROGRAM}
	./${PROGRAM}

# tar the current directory
MY_DIR = `basename ${PWD}`
tar: clean_all
	@echo "Tar the directory: " ${MY_DIR}
	@cd .. ;\
	tar cf ${MY_DIR}.tar ${MY_DIR} *default.mk ;\
	cd ${MY_DIR}
# 	tar cf `basename ${PWD}`.tar *

doc:
	doxygen Doxyfile

#########################################################################

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

##################################################################################################
#    some tools
# Cache behaviour (CXXFLAGS += -g  tracks down to source lines; no -pg in linkflags)
cache: ${PROGRAM}
	valgrind --tool=callgrind --simulate-cache=yes ./$^
#	kcachegrind callgrind.out.<pid> &
	kcachegrind `ls -1tr  callgrind.out.* |tail -1`

# Check for wrong memory accesses, memory leaks, ...
# use smaller data sets
mem: ${PROGRAM}
	valgrind -v --leak-check=yes --tool=memcheck --undef-value-errors=yes --track-origins=yes --log-file=$^.addr.out --show-reachable=yes ./$^

#  Simple run time profiling of your code
#  CXXFLAGS += -g -pg
#  LINKFLAGS += -pg
prof: ${PROGRAM}
	./$^
	gprof -b ./$^ > gp.out
#	kprof -f gp.out -p gprof &
