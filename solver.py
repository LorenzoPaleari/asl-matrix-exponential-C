import sys
import numpy as np
from scipy.linalg import expm
import os

def main():
    # Checks if the script is called with the required parameters
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <tmp_file>")
        return 1
    tmp_file = sys.argv[1]
    
    # Reads the array from the input file
    arr = np.fromfile(tmp_file, dtype=np.float64)
    
    #  delete input_file in the current directory
    os.remove(tmp_file)

    # Reshap the array into a squared matrix
    m = int(np.sqrt(len(arr)))
    matrix = np.reshape(arr, (m, m))

    # Perform the matrix exponential operation
    matrix = expm(matrix)
    
    # Reshapes the matrix into a flattened array
    arr_out = matrix.flatten()

    # Writes the array to the output file
    arr_out.tofile(tmp_file)


if __name__ == "__main__":
    main()