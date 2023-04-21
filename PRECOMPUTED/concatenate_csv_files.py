import os
import pandas as pd

# Define the range of x and y values
x_values = range(101)
y_values = range(101)

# Sort the x and y values
x_values = sorted(x_values)
y_values = sorted(y_values)

# Define the base name of the CSV files
input_base_name = "collated_u{}_d{}_p8_v4.csv"
output_base_name = "collated_p8_v4.csv"

# Define the custom column names
p_value = 8
column_names = ['p_up', 'p_down', 'path'] + [str(i) for i in range(2**p_value + 1)]

# Construct the full path to the input and output files
input_filenames = [os.path.join("MATRICES", input_base_name.format(x, y)) for x in x_values for y in y_values if x + y <= 100]
output_filename = os.path.join("MATRICES", output_base_name)

# Concatenate all CSV files into a single DataFrame
dfs = [pd.read_csv(f, header=None, names=column_names) for f in input_filenames]
concatenated_df = pd.concat(dfs)

# Write the concatenated DataFrame to a new CSV file without a header
concatenated_df.to_csv(output_filename, index=False, header=False)

