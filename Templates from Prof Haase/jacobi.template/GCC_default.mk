# Basic Defintions for using GNU-compiler suite sequentially
# requires setting of COMPILER=GCC_

#startmake as follows to avoid warnings caused by OpenMPI code
#  make 2>&1 | grep -v openmpi


MPI_ROOT=/usr/bin/

CC	= ${MPI_ROOT}mpicc
CXX     = ${MPI_ROOT}mpicxx
F77	= ${MPI_ROOT}mpif77
LINKER  = ${CXX}

MPIRUN  = ${MPI_ROOT}mpirun
#MPIRUN  = ${MPI_ROOT}mpiexec


WARNINGS = -Wall -pedantic -Woverloaded-virtual -Wfloat-equal -Wshadow \
           -Wredundant-decls -Wunreachable-code -Winline -fmax-errors=1

# WARNINGS += -Weffc++ -Wextra
# -Wno-pragmas
CXXFLAGS += -std=c++17 -ffast-math -O3 -march=native ${WARNINGS}
# -ftree-vectorizer-verbose=5  -DNDEBUG
#          -ftree-vectorizer-verbose=2
# CFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp -fdump-tree-vect-details
# CFLAGS	= -ffast-math -O3 -funroll-loops -DNDEBUG -msse3 -fopenmp -ftree-vectorizer-verbose=2

# info on vectorization
#VECTORIZE = -ftree-vectorize -fdump-tree-vect-blocks=foo.dump
#-fdump-tree-pre=stderr
VECTORIZE = -ftree-vectorize -fopt-info -ftree-vectorizer-verbose=5
#CXXFLAGS += ${VECTORIZE}

# -funroll-all-loops   -msse3
#GCC  -march=knl -march=broadwell -march=haswell

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
CXXFLAGS += -fopenmp
LINKFLAGS += -fopenmp

default: ${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@
	@echo
	@echo "Start with :  $(MPIRUN) -np num_proc $(MPIFLAGS) $(PROGRAM)"
	@echo

clean:
	@rm -f ${PROGRAM} ${OBJECTS} gmon.out

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar *.orig
	@rm -rf html latex

run: ${PROGRAM}
	${MPIRUN} -np 4 ./$^

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
#	 2>&1 | grep -v openmpi

# special: get rid of compiler warnings genereate by openmpi-files
.cpp.o:
# 	@$(CXX) -c $(CXXFLAGS) $< 2>/tmp/t.txt || grep -sv openmpi /tmp/t.txt
# 	|grep -sv openmpi

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
	 ${MPIRUN} -np 4 valgrind -v --leak-check=yes --tool=memcheck --undef-value-errors=yes --track-origins=yes --log-file=$^.addr.out --show-reachable=yes ./$^

#  Simple run time profiling of your code
#  CXXFLAGS += -g -pg
#  LINKFLAGS += -pg
prof: ${PROGRAM}
	./$^
	gprof -b ./$^ > gp.out
#	kprof -f gp.out -p gprof &
