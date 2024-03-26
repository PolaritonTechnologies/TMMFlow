import ctypes
import numpy as np
import time
import pandas as pd
import json

import matplotlib.pylab as plt

from scipy.optimize import dual_annealing, minimize, differential_evolution

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

file_name = "calculation_order.json"

my_filter = lib.createFilterStack(file_name.encode("utf-8"))

wavelength = 400
theta_0 = 10.0
phi_0 = 10.0

with open("calculation_order.json") as f:
    data_order = json.load(f)

initial_thicknesses = np.array(data_order["structure_thicknesses"])[::-1]


def merit_function(
    features,
    target_type,
    target_polarization,
    target_value,
    target_condition,
    target_angle,
    target_wavelength,
    target_tolerance,
):
    """
    Function that adds up all the computed values and computes a merit from the
    results
    """
    merit = 0

    for i in range(0, np.size(target_value)):

        # print("my_filter: ", my_filter)
        # print("target_wavelength[i]: ", target_wavelength[i])
        # print("target_angle[i]: ", target_angle[i])
        # print("phi_0: ", phi_0)
        # print("thicknesses: ", thicknesses)
        # print("np.size(thicknesses): ", np.size(thicknesses))

        thicknesses = features.astype(np.float64)
        material_list = np.array([0,1,2], dtype = np.int32)
        
        #lib.change_material_order(my_filter, material_list, np.size(material_list))
        lib.change_material_thickness(my_filter, thicknesses, np.size(thicknesses))
        target_calculated = lib.calculate_reflection_transmission_absorption(
            my_filter,
            target_type[i].encode("utf-8"),
            target_polarization[i].encode("utf-8"),
            float(target_wavelength[i]),
            float(target_angle[i]),
            phi_0,
        )

        if target_condition[i] == "=" and target_calculated != float(target_value[i]):

            merit += (
                (target_calculated - float(target_value[i]))
                / float(target_tolerance[i])
            ) ** 2

        if target_condition[i] == ">" and target_calculated < float(target_value[i]):

            merit += (
                (target_calculated - float(target_value[i]))
                / float(target_tolerance[i])
            ) ** 2

        if target_condition[i] == "<" and target_calculated > float(target_value[i]):

            merit += (
                (target_calculated - float(target_value[i]))
                / float(target_tolerance[i])
            ) ** 2
    print("merit: ", merit)
    print("thicknesses: ", thicknesses)
    return merit


def optimization_function(features):
    """
    Simple optimization function to optimize for a single target
    """
    thicknesses = features[:np.size(features)//2].astype(np.float64)
    material_list = features[np.size(features)//2:].astype(np.int32)
    
    lib.change_material_order(my_filter, material_list, np.size(material_list))
    lib.change_material_thickness(my_filter, thicknesses, np.size(thicknesses))
    reflectivity = lib.calculate_reflection_transmission_absorption(
        my_filter, "r".encode("utf-8"), "s".encode("utf-8"), wavelength, theta_0, phi_0
    )
    return reflectivity


# Test set for position the reflection dip at 450 nm
# target_wavelength = np.array([420,450])
# target = np.array([0.8, 0.6, 0.8])
# tolerance = np.array([1e-3, 1e-3, 1e-3])

with open("optimisation_order.json") as f:
    data_optim = json.load(f)

target_type = np.array(data_optim["targets_type"])
target_polarization = np.array(data_optim["targets_polarization"])
target_value = np.array(data_optim["targets_value"])
target_condition = np.array(data_optim["targets_condition"])
target_angle = np.array(data_optim["targets_angle"])
target_wavelength = np.array(data_optim["targets_wavelengths"])
target_tolerance = np.array(data_optim["targets_tolerance"])
bounds = np.array(data_optim["bounds"])

# a = optimization_function(np.append(initial_thicknesses, [0,1,2]) )
# print(a)

# ret = dual_annealing(merit_function, args = (target, target_wavelength, tolerance), bounds=[(0, 1000), (0, 1000), (0, 1000)])
ret = differential_evolution(
    merit_function,
    args=(
        target_type,
        target_polarization,
        target_value,
        target_condition,
        target_angle,
        target_wavelength,
        target_tolerance,
    ),
    bounds=[(0, 1000), (0, 1000), (0, 1000)],
)
# ret = minimize(merit_function, x0 = initial_thicknesses, args = (target, target_wavelength, tolerance), bounds=[(0, 1000), (0, 1000), (0, 1000)], method = "Nelder-Mead")

# print(ret.x[::-1])

if False:

    initial_time = time.time()
    angles = np.linspace(1, 89, 89)

    # temporary fix for 0 which creates a null kz
    if angles[0] == 0:
        angles[0] = 0.1

    wavelength = np.linspace(400, 700, 301)
    stored_reflectivity = pd.DataFrame(columns=angles.astype("U"), index=wavelength)

    for theta in angles:
        print("Angle: ", theta)
        for wavelength in stored_reflectivity.index:
            stored_reflectivity.loc[
                stored_reflectivity.index == wavelength, str(theta)
            ] = lib.calculate_reflection(
                my_filter,
                wavelength,
                theta,
                phi_0,
                np.ascontiguousarray(initial_thicknesses),
                np.size(initial_thicknesses),
            )

    print(stored_reflectivity)
    print(time.time() - initial_time)
    plt.close()
    X, Y = np.meshgrid(
        stored_reflectivity.columns.astype(float),
        stored_reflectivity.index.astype(float),
    )
    # Prepare data for 3D plot where each column contains the same data for the different angles
    Z = stored_reflectivity.to_numpy(float)
    plt.pcolormesh(X, Y, Z, shading="auto")
    # Save X, Y, Z to csv files
    np.savetxt("X.csv", X, delimiter=",")
    np.savetxt("Y.csv", Y, delimiter=",")
    np.savetxt("Z.csv", Z, delimiter=",")

    # Add a colorbar
    plt.colorbar(label="Reflectivity (%)")

    # Visuals
    plt.xlabel("Angle (°)")
    plt.ylabel("Wavelength (nm)")

    # Save the figure before showing it
    plt.savefig("my_plot.png", format="png", dpi=300)
    plt.show()

    lib.destroyFilterStack(my_filter)
