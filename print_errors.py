import os
from collections import deque

folder_path = "./"
max_line_length = 100
num_last_lines = 5
encoding = "utf-8"

def clean_line(line):
    return "".join(c for c in line if c.isprintable()).rstrip()

for file_name in os.listdir(folder_path):
    if file_name.endswith(".out"):
        file_path = os.path.join(folder_path, file_name)
        print("File:", file_name)
        with open(file_path, "r", encoding=encoding, errors="replace") as file:
            last_lines = deque((clean_line(line) for line in file), maxlen=num_last_lines)
            for line in reversed(last_lines):
                if len(line) <= max_line_length:
                    if "Error" in line:
                        print(line + " [{}]".format(file_name))
                #else:
                #    print("Skipping line:", line)
            #for line in last_lines:
            #    print(line + " [{}]".format(file_name))
