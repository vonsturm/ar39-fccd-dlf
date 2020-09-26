# Makefile
#
# Katharina von Sturm
# 25. September 2020
#
# needs root compiled with -std=c++17

CXX = g++ -Wall -O3
CXXFLAGS = $(shell root-config --cflags) -I.
LIBS = $(shell root-config --libs) -lMinuit -ltbb

EXE = ar39stat

all: $(EXE)

ar39stat: ar39stat.cxx args_reader.cxx arbitrary_sampler.cxx
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@

clean:
	rm -rf $(EXE)