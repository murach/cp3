# CXXFLAGS=-Wall -O3 -openmp -openmp-report1 -g
# LDLIBS=-lgomp
CXXFLAGS=-openmp -openmp-report1
CC=icpc
#icpc

all:	main

main: 	main.cpp frng.cpp
	$(CC) $(CXXFLAGS) -o main main.cpp frng.cpp

# Ndim lambda kappa h
run:	main
	./main 2 2 0.4 0
	./main 2 2 0.55 0

clean:
	rm -f main main.o

.PHONY:	run clean