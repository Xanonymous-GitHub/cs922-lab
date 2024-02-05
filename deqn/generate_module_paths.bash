#!/usr/bin/env bash

# We use this file to put all possible include paths for the compiler
# These include path may be different on different machines.
# Basically, this is a workaround for CLion's inability to parse include paths in a Makefile project.

# Function: get_module_paths
# Description:
#   This function takes expected sub-dir name a list of executable names and generates the corresponding
#   module sub-dir paths. The paths are derived based on the location of the executables.
#   If the derived path is valid and exists, it is printed to stdout.
#
# Arguments:
#   $1 - expected sub-directory name (e.g., include, lib, bin).
#   $@ - List of executable names (e.g., mpic++, g++, g++-13).
#
# Output:
#   Absolute paths to the module sub-directories, one per line.
#
get_module_paths() {
    # Check if the expected sub-directory name is provided
    sub_dir_name="$1"
    if [ -z "$sub_dir_name" ]; then
        # If not, directly return.
        return 0
    fi

    # Loop through each argument provided to the function
    for executable_name in "$@"; do
        # Check if the executable exists in the PATH
        if ! type "$executable_name" > /dev/null 2>&1; then
            continue
        fi

        # Initialize the target link to the executable
        target_link="$(which "$executable_name")"

        # Try to resolve the target link using readlink -f
        resolved_link=$(readlink -f "$target_link" 2>/dev/null)

        if [ -n "$resolved_link" ]; then
            target_link="$resolved_link"
        fi

        # Extract the directory path from the executable path
        dir_path="$(dirname "$target_link")"

        # Use 'realpath' to convert relative paths (../../include) into absolute paths.
        # The '2>/dev/null' part suppresses error messages from 'realpath' in case of failure.
        true_path="$(realpath "$dir_path"/../"$sub_dir_name" 2>/dev/null)"

        # Check if 'true_path' is non-empty and is a directory.
        # If so, print it; otherwise, continue to the next iteration.
        if [ -n "$true_path" ] && [ -d "$true_path" ]; then
            echo "$true_path"
        fi
    done
}

get_module_paths "$1" g++ g++-13
