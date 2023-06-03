CXX = clang++
CXXFLAGS = -O2 -std=c++2a -g 

LIBS = -fopenmp

SRCS = ../gx_random.cpp ../gx_vector.cpp $(wildcard *.cpp)
DEPS = ../gx_random.h ../gx_vector.h ../binary_tree.h $(wildcard *.h)

all: $(SRCS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(WARNINGS) $(SRCS) $(LIBS) -o run
