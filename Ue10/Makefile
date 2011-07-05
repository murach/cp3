# CXXFLAGS=-Wall -O3 -openmp -openmp-report1 -g
# LDLIBS=-lgomp
CXXFLAGS=-O3 -openmp -openmp-report1
CC=icpc
#icpc

all:	main

main: 	main.o frng.o
	$(CC) $(CXXFLAGS) -o main main.cpp frng.cpp

run:	main
	./main

clean:
	rm -f main main.o

.PHONY:	run clean