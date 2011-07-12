# CXXFLAGS=-Wall -O3 -openmp -openmp-report1 -g
# LDLIBS=-lgomp
CXXFLAGS=-openmp -openmp-report1
CC=icpc
#icpc

all:	main

main: 	main.cpp frng.cpp
	$(CC) $(CXXFLAGS) -o main main.cpp frng.cpp

# Ndim lambda kappa h Lmax
# 11.1
run1:	main
	./main 2 2 0.4 0 28
	gnuplot plot.gp
	./main 2 2 0.55 0 128
	gnuplot plot2.gp

# 11.2
run2:	main
	./main 3 2 0.23 0 16
	gnuplot plot3.gp
	./main 3 2 0.3 0 16
	gnuplot plot4.gp

clean:
	rm -f main main.o

.PHONY:	run clean