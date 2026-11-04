# Path to the PROFIL/BIAS installation (override: make PROFIL_DIR=/path/to/Profil-2.0.8)
PROFIL_DIR ?= $(HOME)/programs/Profil-2.0.8

CXX       := g++
CXXFLAGS  := -O0 -w -D__I386__ -I$(PROFIL_DIR)/include -I$(PROFIL_DIR)/include/BIAS
LDFLAGS   := -L$(PROFIL_DIR)/lib
LOADLIBES := -lProfil -lBias -llr -lm

.PHONY: all clean
all: ProblemGenerator-NestedLogit

ProblemGenerator-NestedLogit: ProblemGenerator-NestedLogit.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LOADLIBES)

clean:
	$(RM) ProblemGenerator-NestedLogit
