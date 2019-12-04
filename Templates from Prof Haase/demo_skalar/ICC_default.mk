# Basic Defintions for using INTEL compiler suite sequentially
# requires setting of COMPILER=ICC_

# special on my sony [GH]
#BINDIR = /opt/save.intel/bin/
# very special on my sony [GH]
# FIND_LIBS = -L /opt/save.intel/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_intel_lp64.so

#export KMP_AFFINITY=verbose,compact

CC	= ${BINDIR}icc
CXX     = ${BINDIR}icpc
F77	= ${BINDIR}ifort
LINKER  = ${CXX}

WARNINGS = -pedantic -Wall -Weffc++ -Woverloaded-virtual -Wfloat-equal -Wshadow -wd2015,2012
          #-Winline -Wunreachable-code  -Wredundant-decls
CXXFLAGS +=  -std=c++17 -O3  -fma -DNDEBUG ${WARNINGS} -mkl
#CXXFLAGS +=  -std=c++17 -O3 -march=core-avx2  -fma -ftz -fomit-frame-pointer -DNDEBUG ${WARNINGS} -mkl
# -fast       # fast inludes also -ipo !
CXXFLAGS +=  -fargument-noalias -fargument-noalias-global -ansi-alias
CXXFLAGS +=  -align -qopt-dynamic-align
#CXXFLAGS +=  -xCore-AVX2
#CXXFLAGS +=  -tp=zen
# -qopt-subscript-in-range
# -vec-threshold0
# -xCORE-AVX2
# -axcode COMMON-AVX512 -axcode MIC-AVX512 -axcode CORE-AVX512 -axcode CORE-AVX2
# -ipo

# Reports: https://software.intel.com/en-us/articles/getting-the-most-out-of-your-intel-compiler-with-the-new-optimization-reports
#CXXFLAGS +=  -qopt-report=5 -qopt-report-phase=vec,par

#CXXFLAGS +=  -qopt-report=5 -qopt-report-phase=cg
# Redirect report from *.optrpt to stderr
#    -qopt-report-file=stderr
# Guided paralellization
#    -guide -parallel
#    -guide-opts=string  -guide-par[=n]  -guide-vec[=n]
#    -auto-p32 -simd

# interprocedural optimization
CXXFLAGS += -ipo
LINKFLAGS += -ipo

# OpenMP
CXXFLAGS += -qopenmp
# -qopt-report-phase=openmp
# -diag-enable=sc-full  -diag-file=filename -diag-file-append[=filename]
LINKFLAGS += -qopenmp

# use MKL by INTEL
# LINKFLAGS += -L${BINDIR}../composer_xe_2013.1.117/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
LINKFLAGS += -O2 -mkl
# -ipo



default:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@

clean:
	rm -f ${PROGRAM} ${OBJECTS} *.optrpt

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
	echo 0 | sudo tee /proc/sys/kernel/perf_event_paranoid
	amplxe-gui &

# GUI for Memory and Thread analyzer (race condition)
inspector: ${PROGRAM}
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
	#${BINDIR}../inspector_xe_2013/bin64/inspxe-gui &
	inspxe-gui &

advisor:
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
#	https://software.intel.com/en-us/articles/intel-advisor-2017-update-1-what-s-new
	export ADVIXE_EXPERIMENTAL=roofline
	advixe-gui &

icc-info:
	icpc -# main.cpp




