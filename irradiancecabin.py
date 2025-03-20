import pandas as pd
import pyarrow.parquet as pq
from scipy.io import savemat

# File paths
parquet_file = r"C:\Users\susan\OneDrive - UPV\Desktop_UPV\Doctorado\PhD Simulation Model\Lumped model\2024_Matlab Cabin HVAC model\CabinModel14Nodes\prius_VIPV.parquet"
matlab_file = "prius_VIPV.mat"

# Set chunk size
chunk_size = 1000  # Number of rows to process at a time

# Open the Parquet file with PyArrow
parquet_reader = pq.ParquetFile(parquet_file)

# Dictionary to store data that will be saved to .mat file
mat_data = {}

# Get the total number of row groups in the Parquet file
num_row_groups = parquet_reader.num_row_groups

# Process the file in row groups using pyarrow
try:
    for i in range(num_row_groups):
        # Read the row group and convert to pandas DataFrame
        chunk = parquet_reader.read_row_group(i).to_pandas()

        # If chunk is not empty, add it to the dictionary
        if not chunk.empty:
            mat_data[f"chunk_{i+1}"] = chunk

    # Save the data to a .mat file
    savemat(matlab_file, mat_data)

    print(f"Data has been successfully saved to {matlab_file}")

except Exception as e:
    print(f"An error occurred: {e}")