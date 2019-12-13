# Basic Defintions for using PGI-compiler suite sequentially
# requires setting of COMPILER=PGI_
# OPTIRUN = optirun
# on mephisto:
#CXXFLAGS  += -I/share/apps/atlas/include
#LINKFLAGS += -L/share/apps/atlas/lib
#LINKFLAGS   += -lcblas -latlas

LINKFLAGS   += -lblas

CC	= pgcc
CXX     = pgc++
F77	= pgfortran
LINKER  = ${CXX}


WARNINGS = -Minform=warn

#PGI_PROFILING = -Minfo=loop,vect,opt,intensity,mp,accel
PGI_PROFILING = -Minfo=ccff,accel,ipa,loop,lre,mp,opt,par,unified,vect,intensity

# -Minfo
# -Mprof=lines

CXXFLAGS += -std=c++14 -O3 -fast  -DNDEBUG ${PGI_PROFILING} ${WARNINGS}
CXXFLAGS += -Mvect -Mcache_align -Msafeptr -Mprefetch -Mlre -Mdepchk
#-Msmart  

LINKFLAGS   += ${PGI_PROFILING}
#-lcblas
# OpenMP
CXXFLAGS += -mp=align,bind,numa -Mneginfo=mp
LINKFLAGS += -mp=allcores,bind,numa

default:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@

clean:
	@rm -f ${PROGRAM} ${OBJECTS}

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
# #  Simple run time profiling of your code
# #  CXXFLAGS += -g
# #  LINKFLAGS +=


# Profiling options PGI, see: pgprof -h
PROF_FILE = jac.pgprof
# CPU_PROF = -allcache
CPU_PROF = --cpu-profiling on --analysis-metrics 
# GPU_PROF = -cuda=gmem,branch,cc13 -cudainit
#GPU_PROF = -cuda=branch:cc20
#

cache: prof

prof: ${PROGRAM}
#	./$^
#	$(CUDA_HOME)/bin/nvvp &
#  more  /opt/pgi/linux86-64/16.10/bin/pgcollectrc
	${OPTIRUN} ${BINDIR}pgprof ${CPU_PROF} -o $(PROF_FILE) ./$^
	${OPTIRUN} ${BINDIR}pgprof -i  $(PROF_FILE) 2> prof.out
	

