import ctypes
import json
import numpy as np

from optim_module import OptimModule
from plotting_module import PlottingModule

#########################
# Input parameters
optimisation_order_file = "optimisation.json"

#########################

#########################
# Link C++ functions and create filter
lib = ctypes.CDLL("./run_filter_stack.so")

FilterStack = ctypes.POINTER(ctypes.c_char)

c_double_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

c_int_array = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")

lib.createFilterStack.argtypes = [ctypes.c_char_p]
lib.createFilterStack.restype = FilterStack

lib.destroyFilterStack.argtypes = [FilterStack]

lib.calculate_reflection_transmission_absorption.argtypes = [
    FilterStack,
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
]
lib.calculate_reflection_transmission_absorption.restype = ctypes.c_double

lib.change_material_thickness.argtypes = [FilterStack, c_double_array, ctypes.c_size_t]
lib.change_material_order.argtypes = [FilterStack, c_int_array, ctypes.c_size_t]
lib.reset_filter.argtypes = [FilterStack]

my_filter = lib.createFilterStack(optimisation_order_file.encode("utf-8"))
#########################

#########################
# Optimization
print("running optimisation...")
optimization = OptimModule(optimisation_order_file, my_filter, lib)
features = optimization.perform_optimisation("minimize")
print("optimised features: ", features)
#########################

#########################
# Plotting
plotting = PlottingModule(my_filter, lib)

with open(optimisation_order_file) as f:
    plot_order = json.load(f)

wavelength = np.arange(
    plot_order["wavelengthMin"],
    plot_order["wavelengthMax"] + 1,
    plot_order["wavelengthStep"],
)
polar_angles = np.arange(
    plot_order["polarAngleMin"],
    plot_order["polarAngleMax"] + 1,
    plot_order["polarAngleStep"],
)
azimuthal_angles = np.arange(
    plot_order["azimAngleMin"],
    plot_order["azimAngleMax"] + 1,
    plot_order["azimAngleStep"],
)

print("plotting results...")
plotting.plot_ar_data(wavelength, polar_angles, azimuthal_angles, "r", "s", True)
#########################
