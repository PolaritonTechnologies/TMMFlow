import ctypes
import numpy as np
import time
import pandas as pd
import json
import matplotlib.pylab as plt
import time
from scipy.optimize import dual_annealing, minimize, differential_evolution


class CalcModule:
    def __init__(self, calculation_order_file):

        self.lib = ctypes.CDLL("./run_filter_stack.so")
        FilterStack = ctypes.POINTER(ctypes.c_char)

        c_double_array = np.ctypeslib.ndpointer(
            dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"
        )

        self.lib.createFilterStack.argtypes = [ctypes.c_char_p]
        self.lib.createFilterStack.restype = FilterStack

        self.lib.destroyFilterStack.argtypes = [FilterStack]

        self.lib.calculate_reflection_transmission_absorption.argtypes = [
            FilterStack,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            c_double_array,
            ctypes.c_size_t,
        ]
        self.lib.calculate_reflection_transmission_absorption.restype = ctypes.c_double

        self.my_filter = self.lib.createFilterStack(
            calculation_order_file.encode("utf-8")
        )

        with open(calculation_order_file) as f:
            data_order = json.load(f)

        self.calculation_type = data_order["calculation_type"]
        self.angles = np.arange(
            float(data_order["angleMin"]),
            float(data_order["angleMax"]),
            float(data_order["angleStep"]),
        )
        self.wavelengths = np.arange(
            float(data_order["wavelengthMin"]),
            float(data_order["wavelengthMax"]),
            float(data_order["wavelengthStep"]),
        )
        self.structure_materials = data_order["structure_materials"]
        self.structure_thicknesses = data_order["structure_thicknesses"][::-1]
        self.polarization = data_order["polarization"]
        self.azimuthal_angles = data_order["azimuthalAngles"]

    def perform_calculation(self):
        """
        Function that performs the calculation and returns the results
        """

        stored_reflectivity = pd.DataFrame(
            columns=self.angles.astype("U"), index=self.wavelengths
        )

        for theta in self.angles:
            print("angle: ", theta)
            for wavelength in stored_reflectivity.index:
                stored_reflectivity.loc[
                    stored_reflectivity.index == wavelength, str(theta)
                ] = self.lib.calculate_reflection_transmission_absorption(
                    self.my_filter,
                    self.calculation_type.encode("utf-8"),
                    self.polarization.encode("utf-8"),
                    wavelength,
                    theta,
                    0,
                    np.array(self.structure_thicknesses).astype(np.float64),
                    np.size(self.structure_thicknesses),
                )

        print(stored_reflectivity)

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
        plt.savefig("calculation_plot.png", format="png", dpi=300)
        plt.show()


if __name__ == "__main__":

    # Simple test to see if everything is working
    calc_module = CalcModule("test_calculation.json")
    print(calc_module.perform_calculation())
