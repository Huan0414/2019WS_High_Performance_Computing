#
# use GNU-Compiler tools
COMPILER=g++

# use Intel compilers
#COMPILER=ICC_

# use PGI compilers
# COMPILER=PGI_

run: hello
	./a.out
		
hello:
	${COMPILER} -O3 -fopenmp main.cpp mylib.cpp file_io.cpp
