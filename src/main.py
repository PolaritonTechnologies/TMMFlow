import ctypes
import json
import numpy as np

from optim_module import OptimModule
from calculation_module import CalculationModule

#########################
# Input parameters


def translate_order_for_cpp(optimisation_order_file):

    with open(optimisation_order_file, "r") as optimisation_order_file:
        optimisation_order = json.load(optimisation_order_file)

    updated_optimisation_order = optimisation_order.copy()

    # translate the structure

    structure_materials = optimisation_order["structure_materials"]
    structure_thicknesses = optimisation_order["structure_thicknesses"]
    thickness_opt_allowed = optimisation_order["thickness_opt_allowed"]
    layer_switch_allowed = optimisation_order["layer_switch_allowed"]
    bounds = optimisation_order["bounds"]

    updated_structure_materials = []
    updated_structure_thicknesses = []
    updated_thickness_opt_allowed = []
    updated_layer_switch_allowed = []
    updated_bounds = []

    for m_idx, mat in enumerate(structure_materials):

        ## periodic assemblies start with a number

        if mat[0].isdigit():
            parts = mat.split("_")

            nb_layers = int(parts[0])
            layer1_mat = parts[1]
            layer2_mat = parts[2]

            for i in range(nb_layers):
                updated_structure_materials.append(layer1_mat)
                updated_structure_materials.append(layer2_mat)
                updated_structure_thicknesses.append(structure_thicknesses[m_idx][0])
                updated_structure_thicknesses.append(structure_thicknesses[m_idx][1])
                updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
                updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
                updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
                updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
                updated_bounds.append(bounds[m_idx])
                updated_bounds.append(bounds[m_idx])

        else:
            updated_structure_materials.append(mat)
            updated_structure_thicknesses.append(structure_thicknesses[m_idx])
            updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
            updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
            updated_bounds.append(bounds[m_idx])

    updated_optimisation_order["structure_materials"] = updated_structure_materials
    updated_optimisation_order["structure_thicknesses"] = updated_structure_thicknesses
    updated_optimisation_order["thickness_opt_allowed"] = updated_thickness_opt_allowed
    updated_optimisation_order["layer_switch_allowed"] = updated_layer_switch_allowed
    updated_optimisation_order["bounds"] = updated_bounds

    print(updated_optimisation_order)

    with open("temp_cpp_order.json", "w") as f:
        json.dump(updated_optimisation_order, f)

    return "temp_cpp_order.json"


optimisation_order_file_python = "test_optimisation.json"
optimisation_order_file = translate_order_for_cpp(optimisation_order_file_python)

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
#########################

#########################
# Plotting
plotting = CalculationModule(my_filter, lib)

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
plotting.calculate_ar_data(
    wavelength,
    polar_angles,
    azimuthal_angles,
    "r",
    "s",
    save_figure=True,
    save_data=True,
)
#########################
