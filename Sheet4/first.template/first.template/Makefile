#
#  Compile with
#   make 2>&1 | grep -v openmpi
#  to avoid warnings caused by OpenMPI

# use GNU-Compiler tools
COMPILER=GCC_
# alternatively from the shell
# export COMPILER=GCC_
# or, alternatively from the shell
# make COMPILER=GCC_

MAIN = main
SOURCES = ${MAIN}.cpp greetings.cpp
OBJECTS = $(SOURCES:.cpp=.o)

PROGRAM	= ${MAIN}.${COMPILER}

# uncomment the next to lines for debugging and detailed performance analysis
CXXFLAGS += -g
LINKFLAGS +=

include ../${COMPILER}default.mk
