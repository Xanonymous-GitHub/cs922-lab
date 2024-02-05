#!/usr/bin/env bash

# Define paths relative to the script location
SCRIPT_DIR=$(dirname "$0")
BUILD_EXEC="$SCRIPT_DIR/build/deqn.original"
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

# Execute the command
"$BUILD_EXEC" "$TEST_INPUT"
