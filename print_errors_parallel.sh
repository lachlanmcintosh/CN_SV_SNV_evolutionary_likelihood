#!/usr/bin/env python
#SBATCH --job-name=file_processing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

import os
from collections import deque
import multiprocessing as mp

folder_path = "./"
max_line_length = 100
num_last_lines = 5
encoding = "utf-8"

def clean_line(line):
    return "".join(c for c in line if c.isprintable()).rstrip()

def process_file(file_path):
    file_name = os.path.basename(file_path)
    print("File:", file_name)
    with open(file_path, "r", encoding=encoding, errors="replace") as file:
        last_lines = deque((clean_line(line) for line in file), maxlen=num_last_lines)
        error_line = None
        for line in reversed(last_lines):
            if "Error" in line:
                error_line = line
                break
        if error_line is not None:
            print(error_line + " [{}]".format(file_name))
        #else:
        #    for line in last_lines:
        #        print(line + " [{}]".format(file_name))

if __name__ == "__main__":
    pool = mp.Pool()
    file_paths = [os.path.join(folder_path, file_name) for file_name in os.listdir(folder_path) if file_name.endswith(".out")]
    pool.map(process_file, file_paths)
    pool.close()
    pool.join()
