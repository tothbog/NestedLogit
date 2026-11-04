# Path to the PROFIL/BIAS installation (override: make PROFIL_DIR=/path/to/Profil-2.0.8)
PROFIL_DIR ?= $(HOME)/programs/Profil-2.0.8

SRCDIR := src

SEARCHDIRS := -I$(SRCDIR) -I$(PROFIL_DIR)/include -I$(PROFIL_DIR)/include/BIAS
LIBDIRS    := -L$(PROFIL_DIR)/lib

CXX       := g++
CXXFLAGS  := -O0 -g -w -D__I386__ $(SEARCHDIRS)
LDFLAGS   := $(LIBDIRS)
LOADLIBES := -lProfil -lBias -llr -lm

OBJS := $(addprefix $(SRCDIR)/, NestedLogit.o BLTree.o Box.o GORoutines.o HESS_ari.o \
        XI_ari.o I_util.o InputOutput.o IntervalVectorPlusz.o List.o Problems.o \
        ProcessBox.o Sorting.o)

.PHONY: all clean
all: NestedLogit

NestedLogit: $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp $(wildcard $(SRCDIR)/*.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.C $(wildcard $(SRCDIR)/*.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) NestedLogit
