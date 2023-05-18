CXX=g++
CXXFLAGS=-std=c++17 -O3

DEPS=$(wildcard *.h) $(wildcard *.cpp)
SRCS=$(wildcard *.cpp)
OBJS=lbfgs.o

all: $(DEPS) $(SRCS) $(OBJS)
	g++ -g -O3 $(SRCS)  $(OBJS) -fopenmp -o run

lbfgs.o: ./lib_lbfgs/lbfgs.c
	g++ ./lib_lbfgs/lbfgs.c -c -o lbfgs.o

clean:
	rm -f run lbfgs.o