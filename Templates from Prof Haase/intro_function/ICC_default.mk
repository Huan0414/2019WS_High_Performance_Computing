# Basic Defintions for using INTEL compiler suite sequentially
# requires setting of COMPILER=ICC_

#BINDIR = /opt/intel/bin/

CC	= ${BINDIR}icc
CXX     = ${BINDIR}icpc
F77	= ${BINDIR}ifort
LINKER  = ${CXX}


WARNINGS = -Wall -Weffc++ -Woverloaded-virtual -Wfloat-equal -Wshadow -wd2015,2012 -wn3
#    -Winline -Wredundant-decls -Wunreachable-code
CXXFLAGS +=  -O3 -fargument-noalias -std=c++17 -DNDEBUG ${WARNINGS} -mkl 
# profiling tools
#CXXFLAGS  += -pg
#LINKFLAGS += -pg
# -vec-report=3
# -qopt-report=5 -qopt-report-phase=vec
# -guide -parallel
# -guide-opts=string  -guide-par[=n]  -guide-vec[=n]
# -auto-p32 -simd -msse3 

LINKFLAGS += -O3 

# LAPACK, BLAS: use MKL by INTEL
# LINKFLAGS += -L${BINDIR}../composer_xe_2013.1.117/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
CXXFLAGS  += -mkl
LINKFLAGS += -mkl

# interprocedural optimization
CXXFLAGS  += -ipo
LINKFLAGS += -ipo

default:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@

clean:
	rm -f ${PROGRAM} ${OBJECTS}

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar

run: clean ${PROGRAM}
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
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
#	alternatively to the solution abouve:
            #edit file  /etc/sysctl.d/10-ptrace.conf     and set variable   kernel.yama.ptrace_scope   variable to 0 .
	amplxe-gui &

# GUI for Memory and Thread analyzer (race condition)
inspector: ${PROGRAM}
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
	inspxe-gui &

advisor:
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
	advixe-gui &

icc-info:
	icpc -# main.cpp




