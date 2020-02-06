# Basic Defintions for using OpenMPI with Intel compilers
# requires setting of COMPILER=OPENMPI_ICC_

# Pass Intel Compilers to the OpenMPI wrappers
#    see: https://www.open-mpi.org/faq/?category=mpi-apps#override-wrappers-after-v1.0
EXPORT = export OMPI_CXX=icpc; export OMPI_CC=icc; export OMPI_mpifort=ifort

CC	= mpicc
CXX = mpicxx
F77	= mpifort
LINKER  = ${CXX}

MPIRUN  = ${MPI_BIN}mpirun

WARNINGS = -Wall  -Wextra -pedantic -Woverloaded-virtual  -Wfloat-equal -Wshadow
           #  -Weffc++ -Wunreachable-code -Winline
CXXFLAGS += -fast -fargument-noalias  -DNDEBUG -std=c++17 ${WARNINGS}
CFLAGS   += -O3 -fargument-noalias  -DNDEBUG -Wall  -Wextra -pedantic -Wfloat-equal -Wshadow
# -vec-report=3 -mkl
# -guide -parallel
# -guide-opts=string  -guide-par[=n]  -guide-vec[=n]
# -auto-p32 -simd

# use MKL by INTEL
LINKFLAGS += -O3 -mkl ${MPI_LINK_FLAGS}
# ipo: warning #11021: unresolved __GI_memset
# see:  https://software.intel.com/en-us/articles/ipo-warning-11021-unresolved-symbols-referenced-a-dynamic-library
LINKFLAGS += 

default:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	@( ${EXPORT}; $(LINKER)  $^  ${LINKFLAGS} -o $@ )
	@echo
	@echo "Start with :  $(MPIRUN) -np num_proc $(MPIFLAGS) $(PROGRAM)"
	@echo

clean:
	rm -f ${PROGRAM} ${OBJECTS}

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar

run: ${PROGRAM}
	${MPIRUN} -np 4 ./$^ ${PROG_ARGS}

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
	@( ${EXPORT}; $(CXX) -c $(CXXFLAGS) $< )

.c.o:
	@( ${EXPORT}; $(CC) -c $(CFLAGS) $< )

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
