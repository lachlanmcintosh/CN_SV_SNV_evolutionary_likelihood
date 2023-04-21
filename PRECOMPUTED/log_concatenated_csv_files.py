import pandas as pd
import numpy as np
import re
import os

input_file = "MATRICES/collated_p8_v4.csv"
output_file = "MATRICES/collated_p8_v4_logged.pkl"

# Define a function for processing chunks
def process_chunk(chunk):
    # Generate column names
    p_value = int(np.log2(len(chunk.columns) - 3))

    # Extract the p_value from the input file name
    p_value_from_file = int(re.search(r'p(\d+)', input_file).group(1))

    # Check if the p_value matches the one in the input file name
    if p_value != p_value_from_file:
        raise ValueError(f"The p_value in the file name ({p_value_from_file}) does not match the calculated p_value ({p_value}).")

    column_names = ['p_up', 'p_down', 'path'] + [str(i) for i in range(2**p_value + 1)]

    # Assign column names to the chunk
    chunk.columns = column_names

    # Log transform (natural logarithm) columns 3 onwards
    chunk.iloc[:, 3:] = chunk.iloc[:, 3:].apply(np.log)

    return chunk

# Initialize an empty DataFrame to store the processed chunks
processed_data = pd.DataFrame()

# Calculate the total file size
file_size = os.path.getsize(input_file)

# Define the desired number of chunks
num_chunks = 10^6

# Calculate the approximate chunk size based on file size and desired number of chunks
approx_chunksize_bytes = file_size // num_chunks

# Calculate the approximate number of rows per chunk
header_row_size = len(next(open(input_file, 'r')))
total_rows = sum(1 for row in open(input_file, 'r')) - 1  # Subtract 1 for the header row
approx_rows_per_chunk = (approx_chunksize_bytes // header_row_size) + 1

# Read the CSV file in chunks
chunksize = int(approx_rows_per_chunk)
processed_rows = 0

for chunk in pd.read_csv(input_file, header=None, chunksize=chunksize):
    processed_chunk = process_chunk(chunk)
    processed_data = pd.concat([processed_data, processed_chunk], ignore_index=True)

    processed_rows += len(chunk)
    completion_percentage = (processed_rows / total_rows) * 100
    print(f"Progress: {completion_percentage:.2f}%")

# Save the transformed DataFrame as a pickle file
processed_data.to_pickle(output_file)

