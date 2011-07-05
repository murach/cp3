CXXFLAGS=-Wall -O3 -openmp -openmp-report -g
LDLIBS=-lgomp
CC = icpc
#icpc

all:	main

main: 	main.o frng.o

run:	main
	./main

clean:
	rm -f main main.o

.PHONY:	run clean