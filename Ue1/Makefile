CXXFLAGS=`root-config --cflags` -Wall -fopenmp -g
LDLIBS=`root-config --libs` -lgomp
CC = g++

all:	main

cg:	cg.cpp
	$(CC) -c $< -o cg.o

main: 	main.o

run:	main
	./main

clean:
	rm -f main main.o

.PHONY:	run clean