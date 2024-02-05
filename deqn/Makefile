BUILDDIR := build
OUTDIR := out
PRODUCT := deqn

SRCDIR := src

CXXINCLUDES := $(shell bash generate_module_paths.bash 'include' | sed 's/^/-I/')
CXXLIBS := $(shell bash generate_module_paths.bash 'lib' | sed 's/^/-L/')

HDRS := $(wildcard $(SRCDIR)/*.hpp)

CSRCS := $(wildcard $(SRCDIR)/*.cpp)
FSRCS += $(wildcard $(SRCDIR)/fortran/*.f90)

OBJS := $(CSRCS:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)
OBJS += $(FSRCS:$(SRCDIR)/fortran/%.f90=$(BUILDDIR)/%.o)

BINARY := $(BUILDDIR)/$(PRODUCT)

# g++ flags:
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
    CXX := clang++
else
    CXX := g++
endif
CXXFLAGS_DEBUG := -g -DDEBUG
CXXFLAGS_OPT := -O3 -finline-functions -m64 -funroll-loops
CXXFLAGS_OPENMP := -fopenmp
CXXFLAGS_STANDARD := -std=c++2a
CXXFLAGS_OTHER := -Wall -Wextra -Wpedantic

LD := $(CXX)

# select optimized or debug
CXXFLAGS := $(CXXFLAGS_OPT) $(CXXFLAGS_OTHER) $(CPPFLAGS)
# CXXFLAGS := $(CXXFLAGS_OPT) $(CXXFLAGS_DEBUG) $(CPPFLAGS)

# add openmp flags (comment out for serial build)
CXXFLAGS += $(CXXFLAGS_OPENMP)
LDFLAGS += $(CXXFLAGS_OPENMP)

all : $(BINARY)

$(BINARY) : $(OBJS)
	@echo linking $@
	$(maketargetdir)
	$(LD) $(LDFLAGS) $(CXXFLAGS_STANDARD) -o $@ $^

$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp
	@echo compiling $<
	$(maketargetdir)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_STANDARD) $(CXXINCLUDES) -c -o $@ $<

define maketargetdir
	-@mkdir -p $(dir $@) > /dev/null 2>&1
endef

clean :
	rm -f $(BINARY) $(OBJS)
	rm -rf $(BUILDDIR)
	rm -rf $(OUTDIR)
