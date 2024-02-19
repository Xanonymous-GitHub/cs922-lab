.PHONY: all

BUILDDIR := build
OUTDIR := out
SRCDIR := src

CXXINCLUDES := $(shell bash generate_module_paths.bash 'include' | sed 's/^/-I/')
CXXLIBS := $(shell bash generate_module_paths.bash 'lib' | sed 's/^/-L/')

# g++ flags:
CXX := mpic++
CXXFLAGS_DEBUG := -g -DDEBUG
CXXFLAGS_OPT := -O3 -finline-functions -m64 -funroll-loops
CXXFLAGS_OPENMP := -fopenmp
CXXFLAGS_STANDARD := -std=c++2a
CXXFLAGS_OTHER := -Wall -Wextra -Wpedantic

CXXFLAGS := $(CXXFLAGS_OPT) $(CXXFLAGS_OTHER) $(CPPFLAGS)

# add openmp flags (comment out for serial build)
CXXFLAGS += $(CXXFLAGS_OPENMP)
LDFLAGS += $(CXXFLAGS_OPENMP)

HDRS := $(wildcard $(SRCDIR)/*.hpp)
BINARIES := $(BUILDDIR)/bin2ppm $(BUILDDIR)/diffbin $(BUILDDIR)/pingpong $(BUILDDIR)/colcopy $(BUILDDIR)/karman

define maketargetdir
	-@mkdir -p $(dir $@) > /dev/null 2>&1
endef

define linking
	@printf "\033[92m linking \033[0m\n" $@
	$(maketargetdir)
	$(CXX) $(LDFLAGS) $(CXXLIBS) $(CXXFLAGS_STANDARD) -o $@ $^
endef

define compiling
	@printf "\033[96m compiling \033[0m\n" $<
	$(maketargetdir)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_STANDARD) $(CXXINCLUDES) -c -o $@ $<
endef

# Universal targets

all : $(BINARIES)

clean :
	rm -f $(BINARIES) $(OBJS)
	rm -rf $(BUILDDIR)
	rm -rf $(OUTDIR)

# Compiling the objects

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(HDRS)
	$(compiling)

# Linking the binaries

$(BUILDDIR)/karman: $(BUILDDIR)/alloc.o $(BUILDDIR)/boundary.o $(BUILDDIR)/init.o $(BUILDDIR)/karman.o $(BUILDDIR)/simulation.o
	$(linking)

$(BUILDDIR)/bin2ppm: $(BUILDDIR)/bin2ppm.o $(BUILDDIR)/alloc.o
	$(linking)

$(BUILDDIR)/diffbin: $(BUILDDIR)/diffbin.o
	$(linking)

$(BUILDDIR)/pingpong: $(BUILDDIR)/pingpong.o
	$(linking)

$(BUILDDIR)/colcopy: $(BUILDDIR)/colcopy.o $(BUILDDIR)/alloc.o
	$(linking)
