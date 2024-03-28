# tests.py performs basic unit testing for the optimisation module
# command is pytest tests.py

from optim_module import OptimModule
import ctypes
import numpy as np

#########################
# Input parameters
optimisation_order_file = "test_optimisation.json"

#########################

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

my_filter = lib.createFilterStack(optimisation_order_file.encode("utf-8"))

optimization = OptimModule(optimisation_order_file, my_filter, lib)


def cavity_optimisation():

    thicknesses = optimization.perform_optimisation("minimize")
    return thicknesses


def test():
    assert optimization.check_targets(cavity_optimisation())
