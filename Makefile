CXX=g++
CXXFLAGS=-std=c++17 -O3

DEPS=$(wildcard *.h) $(wildcard *.cpp)
SRCS=$(wildcard *.cpp)

all: $(DEPS)
	g++ -O3 $(SRCS) -fopenmp -o run