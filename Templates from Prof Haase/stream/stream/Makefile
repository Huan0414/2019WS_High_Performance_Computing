CC = gcc
CFLAGS = -O3 
DIMENSIONS = -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 

FF = gfortran
FFLAGS = -O3

all: stream_f.exe stream_c.exe flops.exe

stream_f.exe: stream.f mysecond.o
	$(CC) $(CFLAGS) -c mysecond.c
	$(FF) $(FFLAGS) $(DIMENSIONS)  -c stream.f
	$(FF) $(FFLAGS) stream.o mysecond.o -o stream_f.exe

stream_c.exe: stream.c
	$(CC) $(CFLAGS) $(DIMENSIONS) stream.c -o stream_c.exe

clean:
	rm -f *.exe *.o

# an example of a more complex build line for the Intel icc compiler
stream.icc: stream.c
	icc -O3 -xCORE-AVX2 -ffreestanding -qopenmp -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 stream.c -o stream.omp.AVX2.80M.20x.icc

# GH
flops.exe:  
	$(CC) $(CFLAGS) -DUNIX flops.c -o flops.exe

run: clean all
	./stream_c.exe
	./flops.exe

MY_DIR = `basename ${PWD}`
tar: clean
	@cd .. ;\
	tar cf ${MY_DIR}.tar ${MY_DIR} *default.mk ;\
	cd ${MY_DIR}

# HG
