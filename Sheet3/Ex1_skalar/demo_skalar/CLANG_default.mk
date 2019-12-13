# Basic Defintions for using GNU-compiler suite sequentially
# requires setting of COMPILER=CLANG_

# Additional packes for OPenMP with clang:  libomp5   libomp-dev
# see    http://packages.ubuntu.com/de/source/xenial/misc/openmprtl

#CLANGPATH=/usr/lib/llvm-5.0/bin/
CC	    = ${CLANGPATH}clang
# https://llvm.org/docs/CompileCudaWithLLVM.html
# https://llvm.org/docs/NVPTXUsage.html
#CXX     = ${CLANGPATH}clang++ -lomptarget  -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=/opt/pgi/linux86-64/2017/cuda/8.0
CXX     = ${CLANGPATH}clang++
#F77	= gfortran
LINKER  = ${CXX}

# on mephisto:
#CXXFLAGS  += -I/share/apps/atlas/include
#LINKFLAGS += -L/share/apps/atlas/lib
#LINKFLAGS   += -lcblas -latlas
#LINKFLAGS   += -lblas
# Der <cblas.h> Header muss mit extern "C" versehen werden, damit g++ alles findet.

#http://clang.llvm.org/docs/UsersManual.html#options-to-control-error-and-warning-messages
WARNINGS = -Weverything -Wno-c++98-compat -Wno-sign-conversion -Wno-date-time -Wno-shorten-64-to-32 -Wno-padded  -ferror-limit=1
#-fsyntax-only -Wdocumentation -Wconversion -Wshadow -Wfloat-conversion -pedantic
CXXFLAGS += -O3 -std=c++17 ${WARNINGS}          # don't use -Ofast
# -ftrapv

# interprocedural optimization
CXXFLAGS += -flto
LINKFLAGS += -flto

# OpenMP : compiles, links, but runs only with one thread
#           http://stackoverflow.com/questions/33357029/using-openmp-with-clang
#       -fopenmp-use-tl  #file:///usr/share/doc/clang-3.8-doc/html/UsersManual.html#openmp-features
CXXFLAGS += -fopenmp
#-I/usr/lib/gcc/x86_64-linux-gnu/5/include -Wno-reserved-id-macro -Wno-deprecated
LINKFLAGS += -fopenmp
#CXXFLAGS +=  -fopenmp=libgomp -I/usr/lib/gcc/x86_64-linux-gnu/5/include/  -Wno-reserved-id-macro -Wno-deprecated
#LINKFLAGS += -fopenmp=libgomp -L/usr/lib/gcc/x86_64-linux-gnu/5
#-fopenmp=libomp          # http://stackoverflow.com/questions/33400462/omp-h-file-not-found-when-compiling-using-clang

#   very good check
# http://clang.llvm.org/extra/clang-tidy/
#TIDYFLAGS = -checks='modernize*'
#   good check, see:  http://llvm.org/docs/CodingStandards.html#include-style
TIDYFLAGS = -checks=llvm*,-llvm-header-guard -header-filter=.* -enable-check-profile -extra-arg="-std=c++17" -extra-arg="-fopenmp"
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
