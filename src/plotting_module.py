import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import time


class PlottingModule:
    def __init__(self, my_filter, lib):
        self.my_filter = my_filter
        self.lib = lib

    def plot_ar_data(
        self, wavelength, angles, target_type, polarization, save=False
    ):
        initial_time = time.time()

        # temporary fix for 0 which creates a null kz
        if angles[0] == 0:
            angles[0] = 0.1

        stored_reflectivity = pd.DataFrame(columns=angles.astype("U"), index=wavelength)

        for theta in angles:
            print("Angle: ", theta)
            for wavelength in stored_reflectivity.index:
                stored_reflectivity.loc[
                    stored_reflectivity.index == wavelength, str(theta)
                ] = self.lib.calculate_reflection_transmission_absorption(
                    self.my_filter,
                    target_type.encode("utf-8"),
                    polarization.encode("utf-8"),
                    float(wavelength),
                    float(theta),
                    0,
                )

        # Print time elapsed for the generation of the reflectivity matrix
        print(time.time() - initial_time)

        # Plotting
        plt.close()
        X, Y = np.meshgrid(
            stored_reflectivity.columns.astype(float),
            stored_reflectivity.index.astype(float),
        )
        # Prepare data for 3D plot where each column contains the same data for
        # the different angles
        Z = stored_reflectivity.to_numpy(float)
        plt.pcolormesh(X, Y, Z, shading="auto")

        # Add a colorbar
        plt.colorbar(label="Intensity")

        # Visuals
        plt.xlabel("Angle (°)")
        plt.ylabel("Wavelength (nm)")

        # Save the figure before showing it
        if save:
            plt.savefig("my_plot.png", format="png", dpi=300)
        else:
            plt.show()

        # Save X, Y, Z to csv files
        # np.savetxt("X.csv", X, delimiter=",")
        # np.savetxt("Y.csv", Y, delimiter=",")
        # np.savetxt("Z.csv", Z, delimiter=",")
        # lib.destroyFilterStack(my_filter)
