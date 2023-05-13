CXX=g++
CXXFLAGS=-O3
SRCS=../gx_random.cpp ../gx_vector.cpp
DEPS=../gx_random.h ../gx_vector.h

EXECUTABLE=main

all:
	@echo "Makefile: Choose between" 
	@echo "	Option 1: $$>make colorMatching" 
	@echo "	Option 2: $$>make shrinkImage"

colorMatching: colorMatching.cpp $(SRCS) $(DEPS)
	$(CXX) $(CXXFLAGS) colorMatching.cpp $(SRCS) -fopenmp -o colorMatching

shrinkImage: shrinkImage.cpp $(SRCS) $(DEPS)
	$(CXX) $(CXXFLAGS) shrinkImage.cpp $(SRCS) -fopenmp -o shrinkImage
	
