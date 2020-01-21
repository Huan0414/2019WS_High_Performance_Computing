# Use the MPI-wrappers from the PGI compiler suite.
# requires setting of COMPILER=PGI_MPI_
#
# requires
#          sudo apt install librdmacm1



# Details for run time information
# export PGI_ACC_TIME=1
# unset PGI_ACC_TIME
# export PGI_ACC_NOTIFY=1
# export PGI_ACC_NOTIFY=3
# unset PGI_ACC_NOTIFY


PGI_PATH =  /opt/pgi/linux86-64/2019/bin
#ifeq "$(HOSTNAME)" "mephisto.uni-graz.at"
#  # mephisto
#  PGI_PATH =  /share/apps/pgi/linux86-64/2016/bin
#endif


#MPI_ROOT=${PGI_PATH}mpi/mpich/bin/
MPI_ROOT= ${PGI_PATH}/../mpi/openmpi-3.1.3/bin/
MPIRUN  = ${MPI_ROOT}mpirun

CC  = ${MPI_ROOT}mpicc
CXX = ${MPI_ROOT}mpicxx
#F77 = ${MPI_ROOT}mpif77
ifndef LINKER
  LINKER  = ${CC}
endif
LINKER  = ${CXX}

WARNINGS = -Minform=warn

PGI_PROFILING += -Minfo=loop,vect,opt,intensity,mp,accel
#PGI_PROFILING += -Mprof=lines –Minfo=ccff

CXXFLAGS +=  -e3 -std=c++17 -fast ${PGI_PROFILING} ${WARNINGS} -Mnodepchk
CFLAGS   +=  -fast ${PGI_PROFILING} ${WARNINGS} -Mnodepchk
#
#  for OpenACC
# Target architecture (nvidia,host)
TA_ARCH = host
#TA_ARCH = nvidia,host
#TA_ARCH = -ta=nvidia:cc2+,cuda5.5,fastmath
#TA_ARCH = -acc -DNDEBUG -ta=nvidia:cc2+,cuda5.5,fastmath,keepgpu
#TA_ARCH = -acc -DNDEBUG -ta=nvidia:cc2+,fastmath,keepgpu

#,keepgpu
# CFLAGS = -O3 -ta=$(TA_ARCH)
#CFLAGS   += -B -gopt $(TA_ARCH)
#CXXFLAGS += -B -gopt $(TA_ARCH)
#  -Minfo=all

# libcudart.a is needed for direct CUDA calls
#LINKFLAGS  = -gopt $(TA_ARCH) -L${BINDIR}../lib $(PGI_PROFILING)
# -lcudart

default:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(LINKER)  $^  ${LINKFLAGS} -o $@

clean:
	rm -f ${PROGRAM} ${OBJECTS} *.gpu *gprof.out

clean_all:: clean
	@rm -f *_ *~ *.bak *.log *.out *.tar

#run: clean ${PROGRAM}
run: ${PROGRAM}
	${MPIRUN} -np 4 ${OPTIRUN} ./${PROGRAM}

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
# #  CXXFLAGS += -g -pg
# #  LINKFLAGS += -pg


# Profiling options PGI, see: pgcollect -help
CPU_PROF = -allcache
GPU_PROF = -cuda=gmem,branch,cc13 -cudainit
#GPU_PROF = -cuda=branch:cc20
#
PROF_FILE = pgprof.out

prof: ${PROGRAM}
#	./$^
#	$(CUDA_HOME)/bin/nvvp &
#	export LD_LIBRARY_PATH=/state/partition1/apps/pgi/linux86-64/12.9/lib:$LD_LIBRARY_PATH
	${OPTIRUN} ${BINDIR}pgcollect $(GPU_PROF) ./$^
	${OPTIRUN} ${BINDIR}pgprof -exe ./$^  $(PROF_FILE) &


# Memory checker (slooooow!!!):
# see doc at /usr/local/cuda/doc/cuda-memcheck.pdf
# mem: ${PROGRAM}
# 	$(CUDA_HOME)memcheck ./$^
