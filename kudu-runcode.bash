#!/usr/bin/env bash

set -e

# Function to create a job name and execute a command with `srun`
function run_job() {
    # Default job name prefix
    local default_job_prefix="¡job_"
    local executable=""

    # Get the job name from command line arguments
    while [ "$#" -gt 0 ]; do
        case "$1" in
            --job-name)
                job_name="$2"
                shift 2
                ;;
            *)
                executable="$1"
                shift
                ;;
        esac
    done

    # Check if job name is provided and non-empty
    if [ -z "$job_name" ]; then
        # Generate a job name from the current UTC time and hash it
        job_name="${default_job_prefix}$(date -u +%s | sha256sum | cut -c1-8)"
    fi

    # Check if the executable exists
    if [ -z "$executable" ]; then
        echo ">>> Error: Executable name is required."
        return 1
    elif ! [ -f "$executable" ]; then
        echo ">>> Error: The executable '$executable' does not exist."
        return 1
    fi

    # Convert the executable path to an absolute path
    executable=$(realpath "$executable")

    # Create a temporary SBATCH file
    # Note: the annotation of a SBATCH file is `#SBATCH`, a `#` followed by `SBATCH` without any space.
    sbatch_file="submit.sbatch"
    {
        echo "#!/usr/bin/env bash"
        echo "#SBATCH --job-name=$job_name"
        echo "#SBATCH --account=cs402users"
        echo "#SBATCH --partition=cs402-openmp"
        echo "#SBATCH --nodes=1"
        echo "#SBATCH --ntasks-per-node=1"
        echo "#SBATCH --cpus-per-task=6"
        echo "#SBATCH --time=00:01:00"
        echo "module load cs402-mpi > /dev/null"
        echo "srun $executable"
    } > "$sbatch_file"

    # Submit the job
    echo ">>> Submitting job '$job_name' with executable '$executable'"
    sbatch "$sbatch_file"

    # Remove the temporary SBATCH file
    if [ -f "$sbatch_file" ]; then
        rm "$sbatch_file"
        echo ">>> Temporary file '$sbatch_file' removed."
    fi
}

# Call the function with all passed arguments
# Note: you can use `--job-name` to specify a custom job name.
run_job "$@"
