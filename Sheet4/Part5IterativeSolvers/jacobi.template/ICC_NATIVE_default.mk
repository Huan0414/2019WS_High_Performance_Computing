# Basic Defintions for using INTEL-MPI with its compilers
# requires setting of COMPILER=ICC_NATIVE_

# MPI_ROOT should be defined by shell
#   path to  icpc  is contained in $PATH
MPI_BIN = $(shell dirname `which icpc` | sed 's/bin\/intel64/mpi\/intel64\/bin/g')/
MPI_LIB = $(shell echo ${MPI_BIN} | sed 's/bin/lib/g')

#  Intel-MPI wrappers used gcc as default !!
CC	= ${MPI_BINa}mpicc  -cc=icc
CXX = ${MPI_BIN}mpicxx -cxx=icpc
F77	= ${MPI_BIN}mpif77 -f77=ifort
LINKER  = ${CXX}

MPIRUN  = ${MPI_BIN}mpirun

WARNINGS = -Wall  -Wextra -pedantic -Woverloaded-virtual  -Wfloat-equal -Wshadow
           #  -Weffc++ -Wunreachable-code -Winline
CXXFLAGS += -O3 -fargument-noalias  -DNDEBUG -std=c++17 ${WARNINGS} ${MPI_COMPILE_FLAGS}
CFLAGS   += -O3 -fargument-noalias  -DNDEBUG -Wall  -Wextra -pedantic -Wfloat-equal \
            -Wshadow ${MPI_COMPILE_FLAGS}
# -vec-report=3 -mkl
# -guide -parallel
# -guide-opts=string  -guide-par[=n]  -guide-vec[=n]
# -auto-p32 -simd

# use MKL by INTEL
LINKFLAGS += -mkl ${MPI_LINK_FLAGS}

default:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@
	@echo
	@echo "Start with :  $(MPIRUN) -np num_proc $(MPIFLAGS) $(PROGRAM)"
	@echo

clean:
	rm -f ${PROGRAM} ${OBJECTS}

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar

run: ${PROGRAM}
	(export LD_LIBRARY_PATH=${MPI_LIB}:${LD_LIBRARY_PATH} ;${MPIRUN} -np 4 ./$^ ${PROG_ARGS})

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
# #    some tools
# # Cache behaviour (CXXFLAGS += -g  tracks down to source lines)
# cache: ${PROGRAM}
# 	valgrind --tool=callgrind --simulate-cache=yes ./$^
# #	kcachegrind callgrind.out.<pid> &
# 
# # Check for wrong memory accesses, memory leaks, ...
# # use smaller data sets
# mem: ${PROGRAM}
# 	valgrind -v --leak-check=yes --tool=memcheck --undef-value-errors=yes --track-origins=yes --log-file=$^.addr.out --show-reachable=yes ./$^
# 
# #  Simple run time profiling of your code
# #  CXXFLAGS += -g -pg
# #  LINKFLAGS += -pg
# prof: ${PROGRAM}
# 	./$^
# 	gprof -b ./$^ > gp.out
# #	kprof -f gp.out -p gprof &
# 


mem: inspector
prof: amplifier
cache: amplifier

gap_par_report:
	${CXX}  -c -guide -parallel $(SOURCES) 2> gap.txt
	
# GUI for performance report
amplifier: ${PROGRAM}
	${BINDIR}../vtune_amplifier_xe_2013/bin64/amplxe-gui &

# GUI for Memory and Thread analyzer (race condition)
inspector: ${PROGRAM}
# http://askubuntu.com/questions/41629/after-upgrade-gdb-wont-attach-to-process
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
	${BINDIR}../inspector_xe_2013/bin64/inspxe-gui & 
