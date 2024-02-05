#!/usr/bin/env bash

# Define paths relative to the script location
SCRIPT_DIR=$(dirname "$0")
BUILD_EXEC="$SCRIPT_DIR/build/deqn"
TEST_INPUT="$SCRIPT_DIR/test/square.in"

# Check if the build executable exists
if [ ! -f "$BUILD_EXEC" ]; then
    echo ">>>> Error: The build executable does not exist at $BUILD_EXEC."
    exit 1
fi

# Check if the test input file exists
if [ ! -f "$TEST_INPUT" ]; then
    echo ">>>> Error: The test input file does not exist at $TEST_INPUT."
    exit 1
fi

# Dynamically set OMP_NUM_THREADS based on the system's available CPU threads
if [ "$(uname)" = "Linux" ]; then
    NUM_CPUS=$(lscpu | grep -E '^CPU\(s\):' | awk '{print $2}')
    THREADS_PER_CORE=$(lscpu | grep -E '^Thread\(s\) per core:' | awk '{print $4}')
    NUM_THREADS=$((NUM_CPUS * THREADS_PER_CORE))
elif [ "$(uname)" = "Darwin" ]; then
    NUM_THREADS=$(sysctl -a | grep machdep.cpu.thread_count | awk '{print $2}')
else
    echo "checking thread numbers: Unsupported operating system."
    NUM_THREADS=1
fi

# Set OpenMP environment variables
export OMP_NUM_THREADS=$NUM_THREADS
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# Execute the command
"$BUILD_EXEC" "$TEST_INPUT"
