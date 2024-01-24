# Define the compiler
CXX := clang++

# Compiler flags
CXXFLAGS := -O -std=c++2b -fopenmp

# Source files - list all .cpp files in the current directory
SOURCES := $(wildcard *.cpp)

# Build directory
BUILDDIR := build

# Executable files, one for each source file, in the build directory
EXECUTABLES := $(patsubst %.cpp,$(BUILDDIR)/%,$(SOURCES))

# Default target
all: $(BUILDDIR) $(EXECUTABLES)

# Rule to ensure the build directory exists
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Rule to compile each source file to its own executable in the build directory
$(BUILDDIR)/%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Clean target
clean:
	rm -rf $(BUILDDIR)

# PHONY targets
.PHONY: all clean
