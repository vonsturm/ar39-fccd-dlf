# Makefile
#
# Katharina von Sturm
# 25. September 2020
#
# needs root compiled with -std=c++17

CXX = g++ -Wall -O3 -Wno-comments -Wno-unused-result
CXXFLAGS = $(shell root-config --cflags) -I./utils -I.
LIBS = $(shell root-config --libs) -lMinuit -ltbb

BINDIR = ./bin
UTILS = $(wildcard ./utils/*.hpp)

EXE = $(BINDIR)/ar39stat

all: $(BINDIR) $(EXE)

$(BINDIR)/ar39stat: ar39stat.cxx $(UTILS)
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@

$(BINDIR):
	mkdir -p $@

clean:
	rm -rf $(EXE)

dist-clean: clean
	rm -rf $(BINDIR)
			