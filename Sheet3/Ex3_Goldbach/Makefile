#
# use GNU-Compiler tools
COMPILER=g++

# use Intel compilers
#COMPILER=ICC_

# use PGI compilers
# COMPILER=PGI_

run: hello
	./a.out
		
hello: main.cpp mylib.h get_primes.cpp count_goldbach.cpp
#	${COMPILER} -O3 -fopenmp main.cpp mylib.h get_primes.cpp count_goldbach.cpp
	${COMPILER} -O3 -fopenmp  $^
