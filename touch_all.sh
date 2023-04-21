#!/bin/bash
#SBATCH --job-name=touch_all
#SBATCH --output=touch_all.out
#SBATCH --error=touch_all.err
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Load any necessary modules (e.g., module load findutils)

# Change to the target directory
cd /vast/scratch/users/lmcintosh/

# Path to a marker file
marker_file="touch_all.marker"

# Check if the marker file exists and has been modified within the last 24 hours
if [ -e "$marker_file" ] && [ "$(find "$marker_file" -mtime -1)" ]; then
    echo "The script has already been executed within the last 24 hours. Exiting."
    exit 0
fi

# Execute the find command with touch on all files and directories
find . -type f -exec touch {} \;
find . -type d -exec touch {} \;

# Update the marker file's timestamp
touch "$marker_file"

# End of the Slurm script

