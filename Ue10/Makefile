# CXXFLAGS=-Wall -O3 -openmp -openmp-report1 -g
# LDLIBS=-lgomp
CXXFLAGS=-openmp -openmp-report1
CC=icpc
#icpc

all:	main

main: 	main.cpp frng.cpp
	$(CC) $(CXXFLAGS) -o main main.cpp frng.cpp

run:	main
	./main

clean:
	rm -f main main.o

.PHONY:	run clean