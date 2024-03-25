#!/usr/bin/env bash

# Define paths relative to the script location
SCRIPT_DIR=$(dirname "$0")
BUILD_EXEC="$SCRIPT_DIR/build/karman"

# Check if the build executable exists
if [ ! -f "$BUILD_EXEC" ]; then
    echo ">>>> Error: The build executable does not exist at $BUILD_EXEC."
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
export MPIRUN_OPTIONS="-binding "domain=omp:compact" -print-rank-map -envall"
export KMP_AFFINITY=compact,1,0,granularity=fine
export I_MPI_HYDRA_BRANCH_COUNT=-1

# Execute the command
mpirun -n 2 "$BUILD_EXEC" -v 2
