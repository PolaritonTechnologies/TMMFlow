import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import math
import time


def ignore(msg=None):
    pass


class CalculationModule:
    def __init__(
        self,
        my_filter,
        lib,
        plot_queue=None,
        plot_done_queue=None,
        log_plot=ignore,
        log_plot_done=ignore,
        core_selection=None,
        web=True,
    ):
        self.my_filter = my_filter
        self.lib = lib
        self.plot_queue = plot_queue
        self.plot_done_queue = plot_done_queue
        self.log_plot = log_plot
        self.log_plot_done = log_plot_done
        self.web = web
        self.core_selection = core_selection

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

        if self.core_selection == None:
            if polarization != "s":
                is_general_core = self.lib.getGeneralMaterialsInStack(self.my_filter)
            else:
                is_general_core = False
        else:
            if self.core_selection == "general":
                is_general_core = True
            if self.core_selection == "fast":
                is_general_core = False

        initial_time = time.time()
        stored_value = pd.DataFrame(columns=polar_angles.astype("U"), index=wavelength)
        for phi in azim_angles:
            for theta in polar_angles:
                print("Polar angle: ", theta, "Azimuthal angle: ", phi)
                for wavelength in stored_value.index:
                    stored_value.loc[stored_value.index == wavelength, str(theta)] = (
                        self.lib.calculate_reflection_transmission_absorption(
                            self.my_filter,
                            target_type.encode("utf-8"),
                            polarization.encode("utf-8"),
                            wavelength,
                            theta,
                            phi,
                            is_general_core,
                        )
                    )

                self.log_plot(
                    {
                        "intensity": np.nan_to_num(
                            stored_value.to_numpy(float)
                        ).tolist(),
                        "wavelength": stored_value.index.astype(float).tolist(),
                        "angles": stored_value.columns.astype(float).tolist(),
                        "azimuthal_angle": phi,
                    }
                )

        self.log_plot_done()

        if not self.web or True:
            # Print time elapsed for the generation of the reflectivity matrix
            print(time.time() - initial_time)
            # Plotting
            plt.close()
            X, Y = np.meshgrid(
                stored_value.columns.astype(float),
                stored_value.index.astype(float),
            )
            # Prepare data for 3D plot where each column contains the same data for
            # the different angles
            Z = stored_value.to_numpy(float)
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
                with open("value.csv", "w") as the_file:
                    the_file.write("\n".join(header_lines))

                # Save actual data by appending
                stored_value.to_csv("value.csv", sep=",", header=True, mode="a")

            plt.show()
