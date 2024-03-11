import ctypes
import os

# Load the shared library
example = ctypes.CDLL('./mysum.so')

# Define the argument types and return type of the function
example.add.argtypes = [ctypes.c_int, ctypes.c_int]
example.add.restype = ctypes.c_int

# Call the function
result = example.add(2, 3)
print(result)  # Outputs: 5