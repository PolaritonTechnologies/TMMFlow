import ctypes
import numpy as np

# Load the shared library
example = ctypes.CDLL('./src/main.so')

# Define the argument types and return type of the function
example.run_function.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
example.run_function.restype = ctypes.c_double

# Define the parameters
wavelength = 400.0
theta_0 = 0.0017453292519943296
phi_0 = 0.0

# Call the function with the parameters
reflection =  example.run_function(wavelength, theta_0, phi_0)

print(9)

# m_r_ps_out and m_t_ps_out now contain the output matrices