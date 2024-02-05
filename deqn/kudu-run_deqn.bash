#!/usr/bin/env bash

set -e

# Function to create a job name and execute a command with `srun`
function run_job() {
    # Default job name prefix
    local default_job_prefix="¡"
    local executable="./run_deqn_square.bash"
    local executable_original="./run_deqn_original_square.bash"

    job_name="${default_job_prefix}$(date -u +%s | sha256sum | cut -c1-8)"

    # Check if the executable exists
    if [ -z "$executable" ]; then
        echo ">>> Error: Executable name is required."
        return 1
    elif ! [ -f "$executable" ]; then
        echo ">>> Error: The executable '$executable' does not exist."
        return 1
    fi

    # Check if the executable_original exists
    if [ -z "$executable_original" ]; then
        echo ">>> Error: Executable name is required."
        return 1
    elif ! [ -f "$executable_original" ]; then
        echo ">>> Error: The executable '$executable_original' does not exist."
        return 1
    fi

    # Convert the executable path to an absolute path
    executable=$(realpath "$executable")

    # The name of the expect exist original version of executable
    executable_original=$(realpath "$executable_original")

    # Create a temporary SBATCH file
    # Note: the annotation of a SBATCH file is `#SBATCH`, a `#` followed by `SBATCH` without any space.
    sbatch_file="submit.sbatch"
    {
        echo "#!/usr/bin/env bash"
        echo "#SBATCH --job-name=$job_name"
        echo "#SBATCH --account=cs402users"
        echo "#SBATCH --partition=cs402-openmp"
        echo "#SBATCH --nodes=2"
        echo "#SBATCH --ntasks-per-node=1"
        echo "#SBATCH --cpus-per-task=40"
        echo "#SBATCH --time=00:01:00"
        echo "#SBATCH --wait-all-nodes=1"
        echo "module load cs402-mpi > /dev/null"
        echo "module load GCC/12.2.0 > /dev/null"
        echo "srun --ntasks=1 -N1 -r 0 --exclusive $executable &"
        echo "srun --ntasks=1 -N1 -r 1 --exclusive $executable_original &"
        echo "wait"
    } > "$sbatch_file"

    # Submit the job
    echo ">>> Submitting job '$job_name' with executable '$executable' and '$executable_original'"
    sbatch "$sbatch_file"

    # Remove the temporary SBATCH file
    if [ -f "$sbatch_file" ]; then
        rm "$sbatch_file"
        echo ">>> Temporary file '$sbatch_file' removed."
    fi
}

# Call the function with all passed arguments
run_job