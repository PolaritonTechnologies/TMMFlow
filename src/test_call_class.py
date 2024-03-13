import ctypes
import numpy as np
import time
import pandas as pd

import matplotlib.pylab as plt

from scipy.optimize import dual_annealing, minimize, differential_evolution

lib = ctypes.CDLL("./run_filter_stack.so")
FilterStack = ctypes.POINTER(ctypes.c_char)

c_double_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

lib.createFilterStack.argtypes = [ctypes.c_char_p]
lib.createFilterStack.restype = FilterStack

lib.destroyFilterStack.argtypes = [FilterStack]

lib.calculate_reflection.argtypes = [
    FilterStack,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    c_double_array,
    ctypes.c_size_t,
]
lib.calculate_reflection.restype = ctypes.c_double


file_name = "calculation_order.json"

my_filter = lib.createFilterStack(file_name.encode("utf-8"))

wavelength = 400
theta_0 = 0.0017453292519943296
phi_0 = 0.0

initial_thicknesses = np.array([100, 105, 20], dtype=np.float64)


def merit_function(thicknesses, target, target_wavelength, tolerance):
    """
    Function that adds up all the computed values and computes a merit from the
    results
    """
    merit = 0
    for i in range(0, np.size(target)):
        reflectivity = lib.calculate_reflection(
            my_filter,
            target_wavelength[i],
            theta_0,
            phi_0,
            thicknesses,
            np.size(thicknesses),
        )
        merit += ((reflectivity - target[i]) / tolerance[i]) ** 2

    print(merit)
    return merit


def optimization_function(thicknesses):
    """
    Simple optimization function to optimize for a single target
    """
    reflectivity = lib.calculate_reflection(
        my_filter, wavelength, theta_0, phi_0, thicknesses, np.size(thicknesses)
    )
    return reflectivity


target_wavelength = np.array([400, 500, 600])
target = np.array([0.8, 0.6, 0.8])
tolerance = np.array([1e-3, 1e-3, 1e-3])

# ret = dual_annealing(merit_function, args = (target, target_wavelength, tolerance), bounds=[(0, 1000), (0, 1000), (0, 1000)])
# ret = differential_evolution(merit_function, args = (target, target_wavelength, tolerance), bounds=[(0, 1000), (0, 1000), (0, 1000)])
# ret = minimize(merit_function, x0 = initial_thicknesses, args = (target, target_wavelength, tolerance), bounds=[(0, 1000), (0, 1000), (0, 1000)], method = "Nelder-Mead")

# print(ret)

initial_time = time.time()
angles = np.arange(0, 89, 1)
wavelength = np.arange(250, 1000, 1)
stored_reflectivity = pd.DataFrame(columns = angles.astype("U"), index = wavelength)

for theta in angles:
    print("Angle: ", theta)
    for wavelength in stored_reflectivity.index:
        stored_reflectivity.loc[stored_reflectivity.index == wavelength, str(theta)] = lib.calculate_reflection(
            my_filter,
            wavelength,
            theta,
            phi_0,
            initial_thicknesses,
            np.size(initial_thicknesses),
        )

print(stored_reflectivity)
print(time.time() - initial_time)
plt.close()
X, Y = np.meshgrid(stored_reflectivity.columns.astype(float), stored_reflectivity.index.astype(float))
# Prepare data for 3D plot where each column contains the same data for the different angles
Z = stored_reflectivity.to_numpy(float) 
plt.pcolormesh(X, Y,Z , shading="auto")

# Add a colorbar
plt.colorbar(label = "Reflectivity (%)")

# Visuals
plt.xlabel("Angle (°)")
plt.ylabel("Wavelength (nm)")
# plt.ylim(300, 2000)

plt.show()

lib.destroyFilterStack(my_filter)
