#!/usr/bin/env bash

# Build the original version of the code
git checkout original_deqn
git pull
make clean && make

# Build the optimized version of the code
git checkout main
git pull
make
