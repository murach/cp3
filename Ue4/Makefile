CXXFLAGS=`root-config --cflags` -Wall -O3 -fopenmp -g
LDLIBS=`root-config --libs` -lgomp
CC = g++

all:	main

main: 	main.o

run:	main
	./main

clean:
	rm -f main main.o

.PHONY:	run clean