#!/usr/bin/env bash

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

    # Convert the executable path to an absolute path
    if [ -z "$executable" ]; then
        echo ">>> Error: Executable name is required."
        return 1
    elif ! executable=$(realpath "$executable"); then
        echo ">>> Error: The executable '$executable' does not exist or cannot be accessed."
        return 1
    fi

    # Create a temporary SBATCH file
    sbatch_file="submit.sbatch"
    {
        echo "#!/usr/bin/env bash"
        echo "sbatch --job-name=$job_name"
        echo "sbatch --account=cs402users"
        echo "sbatch --partition=cs402-openmp"
        echo "sbatch --nodes=1"
        echo "sbatch --ntasks-per-node=1"
        echo "sbatch --cpus-per-task=6"
        echo "sbatch --time=00:01:00"
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

# Load the MPI module
module load cs402-mpi

# Call the function with all passed arguments
# Note: you can use `--job-name` to specify a custom job name.
run_job "$@"
