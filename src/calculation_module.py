import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import time


class CalculationModule:
    def __init__(self, my_filter, lib):
        self.my_filter = my_filter
        self.lib = lib

    def calculate_ar_data(
        self,
        wavelength,
        polar_angles,
        azim_angles,
        target_type,
        polarization,
        save_figure=False,
        save_data=False,
    ):

        initial_time = time.time()
        stored_reflectivity = pd.DataFrame(
            columns=polar_angles.astype("U"), index=wavelength
        )
        for phi in azim_angles:
            for theta in polar_angles:
                print("Polar angle: ", theta, "Azimuthal angle: ", phi)
                for wavelength in stored_reflectivity.index:
                    stored_reflectivity.loc[
                        stored_reflectivity.index == wavelength, str(theta)
                    ] = self.lib.calculate_reflection_transmission_absorption(
                        self.my_filter,
                        target_type.encode("utf-8"),
                        polarization.encode("utf-8"),
                        wavelength,
                        theta,
                        phi,
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
            plt.colorbar(label="Intensity (a.u.)")
            # Visuals
            plt.xlabel("Polar Angle (°)")
            plt.ylabel("Wavelength (nm)")
            plt.title(f"Calculation for Azimuthal Angle {phi}°")
            # Save the figure before showing it

            if save_figure:
                # plt.savefig(f"{phi}-plot.png", format="png", dpi=300)
                plt.savefig("plot.png", format="png", dpi=300)
                # Save X, Y, Z to csv files
            if save_data:
                header_lines = []
                header_lines.append("Here we can add some meta data to the header \n")

                # Save header lines indicating what the simulation represents
                with open("reflectivity.csv", "w") as the_file:
                    the_file.write("\n".join(header_lines))
                
                # Save actual data by appending
                stored_reflectivity.to_csv("reflectivity.csv", sep = ",", header = True, mode = "a")

                print("Reflectivity data saved!")

            plt.show()
