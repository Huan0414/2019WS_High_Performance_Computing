run: hello
	./hello
		

hello: Sheet2Ex3.o FunctionsSheet2Ex3.o
	g++ Sheet2Ex3.o FunctionsSheet2Ex3.o -o hello

Sheet2Ex3.o: Sheet2Ex3.cc functionsEx3Sheet2.h
	g++ -c Sheet2Ex3.cc

FunctionsSheet2Ex3.o: FunctionsSheet2Ex3.cc functionsEx3Sheet2.h
	g++ -c FunctionsSheet2Ex3.cc

clean:
	rm -rf *o Run_It
