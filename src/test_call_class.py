import ctypes
import numpy as np

lib = ctypes.CDLL('./run_filter_stack.so')
FilterStack = ctypes.POINTER(ctypes.c_char)

c_double_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')

lib.createFilterStack.argtypes = [ctypes.c_char_p]
lib.createFilterStack.restype = FilterStack

lib.destroyFilterStack.argtypes = [FilterStack]

lib.calculate_reflection.argtypes = [FilterStack, ctypes.c_double, ctypes.c_double, ctypes.c_double, c_double_array, ctypes.c_size_t]
lib.calculate_reflection.restype = ctypes.c_double


file_name = "calculation_order.json" 

my_filter = lib.createFilterStack(file_name.encode('utf-8'))

wavelength = 400
theta_0 = 0.0017453292519943296
phi_0 = 0.0

initial_thicknesses = np.array([100, 105, 20], dtype=np.float64)

test = lib.calculate_reflection(my_filter, wavelength, theta_0, phi_0, initial_thicknesses, np.size(initial_thicknesses))
print(test)

lib.destroyFilterStack(my_filter)