#!/usr/bin/env sh

# We use this file to put all possible include paths for the compiler
# These include path may be different on different machines.
# Basically, this is a workaround for CLion's inability to parse include paths in a Makefile project.

# Function: get_include_paths
# Description:
#   This function takes a list of executable names and generates the corresponding
#   include paths. The paths are derived based on the location of the executables.
#   If the derived path is valid and exists, it is printed to stdout.
#
# Arguments:
#   $@ - List of executable names (e.g., mpic++, g++, g++-13).
#
# Output:
#   Absolute paths to the include directories, one per line.
#
get_include_paths() {
    # Loop through each argument provided to the function
    for executable_name in "$@"; do
        # Check if the executable exists in the PATH
        if ! type "$executable_name" > /dev/null 2>&1; then
            continue
        fi

        # Find the absolute path of the executable using 'which',
        # then use 'readlink' to resolve any symbolic links.
        # 'path' will contain the resolved absolute path to the executable.
        path="$(readlink -f "$(which "$executable_name")")"

        # Use 'realpath' to convert relative paths (../../include) into absolute paths.
        # The '2>/dev/null' part suppresses error messages from 'realpath' in case of failure.
        true_path="$(realpath "$path"/../../include 2>/dev/null)"

        # Check if 'true_path' is non-empty and is a directory.
        # If so, print it; otherwise, continue to the next iteration.
        if [ -n "$true_path" ] && [ -d "$true_path" ]; then
            echo "$true_path"
        fi
    done
}

get_include_paths mpic++ g++ g++-13
