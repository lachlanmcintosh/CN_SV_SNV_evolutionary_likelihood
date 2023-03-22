import os

folder_path = "./"
max_line_length = 100
num_last_lines = 5
encoding = "utf-8" # Change this to the correct encoding if needed

def clean_line(line):
    # Remove any non-printable characters and trailing whitespace
    return "".join(c for c in line if c.isprintable()).rstrip()

for file_name in os.listdir(folder_path):
    if file_name.endswith(".out"):
        file_path = os.path.join(folder_path, file_name)
        print("File:", file_name)
        with open(file_path, "r", encoding=encoding, errors="replace") as file:
            lines = file.readlines()
            last_lines = []
            for line in reversed(lines):
                line = clean_line(line)
                if len(line) <= max_line_length:
                    last_lines.append(line)
                    if len(last_lines) == num_last_lines:
                        break
                #else:
                #    print("Skipping line:", line)
            last_lines.reverse()
            print("Last lines:")
            for line in last_lines:
                print(line)
