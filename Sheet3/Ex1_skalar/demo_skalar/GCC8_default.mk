# Basic Defintions for using GNU-compiler suite sequentially
# requires setting of COMPILER=GCC_

CC	= gcc
CXX     = /usr/bin/g++-8
F77	= gfortran
LINKER  = ${CXX}

# on mephisto:
#CXXFLAGS  += -I/share/apps/atlas/include
#LINKFLAGS += -L/share/apps/atlas/lib -L/usr/lib64/atlas
#LINKFLAGS   += -latlas -lcblas
#LINKFLAGS   += -lblas
# Der <cblas.h> Header muss mit extern "C" versehen werden, damit g++ alles findet.


WARNINGS = -Wall -pedantic -Wextra -Weffc++ -Woverloaded-virtual -Wfloat-equal -Wshadow \
           -Wredundant-decls -Winline -fmax-errors=1
#  -Wunreachable-code
#CXXFLAGS += -std=c++17 -ffast-math -O3 -march=native -DNDEBUG ${WARNINGS}
CXXFLAGS += -std=c++17 -ffast-math -O3 -march=native ${WARNINGS}

# info on vectorization
#VECTORIZE = -ftree-vectorize -fdump-tree-vect-blocks=foo.dump
#-fdump-tree-pre=stderr
VECTORIZE = -ftree-vectorize -fopt-info -ftree-vectorizer-verbose=5
#CXXFLAGS += ${VECTORIZE}
# -funroll-all-loops   -msse3
#GCC  -march=knl -march=broadwell -march=haswell

# interprocedural optimization
#CXXFLAGS += -flto
LINKFLAGS += -flto

# for debugging purpose (save code)
# -fsanitize=leak         # only one out the trhee can be used
# -fsanitize=address
# -fsanitize=thread
SANITARY =  -fsanitize=address  -fsanitize=undefined -fsanitize=null -fsanitize=return \
 -fsanitize=bounds -fsanitize=alignment -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow \
 -fsanitize=bool -fsanitize=enum -fsanitize=vptr
#CXXFLAGS  += ${SANITARY}
#LINKFLAGS +=${SANITARY}

# OpenMP
CXXFLAGS += -fopenmp -DOLD_OPENMP
LINKFLAGS += -fopenmp

default: ${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@

clean:
	@rm -f ${PROGRAM} ${OBJECTS}

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar *.orig *.optrpt
	@rm -rf html

run: clean ${PROGRAM}
#	time  ./${PROGRAM}
	./${PROGRAM}

# tar the current directory
MY_DIR = `basename ${PWD}`
tar:
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


thread:${PROGRAM}
	valgrind -v --tool=helgrind --log-file=$^.thread.out ./$^

#  Simple run time profiling of your code
#  CXXFLAGS += -g -pg
#  LINKFLAGS += -pg
prof: ${PROGRAM}
	./$^
	gprof -b ./$^ > gp.out
#	kprof -f gp.out -p gprof &
