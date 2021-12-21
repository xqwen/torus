CXX=g++
CXXFLAGS+=-O3 -Wall
LDLIBS+=-lm -lgsl -lgslcblas -lboost_iostreams -lz

UNAME_S := $(shell uname -s)

# MacOS / Homebrew specific configuration
ifeq ($(UNAME_S),Darwin)
	CPPFLAGS+=-isystem /opt/homebrew/include # Mac libraries installed via Homebrew
	CPPFLAGS+=-Xpreprocessor -fopenmp # Mac openmp
	LDFLAGS+=-L/opt/homebrew/lib # Homebrew library location
	LDLIBS+=-lomp #only need lomp on Mac OS
endif

objects = main.o controller.o logistic.o

all: main

main: $(objects)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(objects) $(LDLIBS) -o torus
static: $(objects)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(objects) $(LDLIBS) -static -o torus.static
main.o: main.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c  main.cc
controller.o: controller.cc classdef.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c  controller.cc
logistic.o: logistic.cc logistic.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c logistic.cc -fpermissive
clean:
	rm *.o torus	
