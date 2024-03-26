# import ctypes
import numpy as np
import json
from scipy.optimize import dual_annealing, minimize, differential_evolution


class OptimModule:
    def __init__(self, optimisation_order_file, my_filter, ctypes_lib):

        self.lib = ctypes_lib
        self.my_filter = my_filter

        with open(optimisation_order_file) as f:
            self.data_order = json.load(f)

        self.initial_thicknesses = np.array(self.data_order["structure_thicknesses"])[::-1]
        self.target_type = np.array(self.data_order["targets_type"])
        self.target_polarization = np.array(self.data_order["targets_polarization"])
        self.target_value = np.array(self.data_order["targets_value"])
        self.target_condition = np.array(self.data_order["targets_condition"])
        self.target_angle = np.array(self.data_order["targets_angle"])
        self.target_wavelength = np.array(self.data_order["targets_wavelengths"])
        self.target_tolerance = np.array(self.data_order["targets_tolerance"])
        self.bounds = [tuple(x) for x in self.data_order["bounds"]]

    def merit_function(self, thicknesses):
        """
        Function that adds up all the computed values and computes a merit from the
        results
        """
        merit = 0

        for i in range(0, np.size(self.target_value)):

            self.lib.change_material_thickness(
                self.my_filter, thicknesses, np.size(thicknesses)
            )
            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                float(self.target_wavelength[i]),
                float(self.target_angle[i]),
                0,
            )

            print(self.target_wavelength[i], target_calculated)

            if self.target_condition[i] == "=" and target_calculated != float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - float(self.target_value[i]))
                    / float(self.target_tolerance[i])
                ) ** 2

            if self.target_condition[i] == ">" and target_calculated < float(
                self.target_value[i]
            ):
                print(
                    self.target_wavelength[i],
                    "> : ",
                    target_calculated,
                    float(self.target_value[i]),
                )

                merit += (
                    (target_calculated - float(self.target_value[i]))
                    / float(self.target_tolerance[i])
                ) ** 2

            if self.target_condition[i] == "<" and target_calculated > float(
                self.target_value[i]
            ):
                print(
                    self.target_wavelength[i],
                    "< : ",
                    target_calculated,
                    float(self.target_value[i]),
                )

                merit += (
                    (target_calculated - float(self.target_value[i]))
                    / float(self.target_tolerance[i])
                ) ** 2

        print("merit: ", merit)
        print("thicknesses: ", thicknesses)
        return merit

    # check targets is used in unit testing
    def check_targets(self, thicknesses):
        """
        Function that returns true if all targets are met, for unit testing only
        """
        test_pass = True

        for i in range(0, np.size(self.target_value)):

            self.lib.change_material_thickness(self.my_filter, thicknesses, np.size(thicknesses))

            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                float(self.target_wavelength[i]),
                float(self.target_angle[i]),
                0,
            )

            if self.target_condition[i] == "=" and target_calculated != float(
                self.target_value[i]
            ):

                print(
                    "Target not met: at ",
                    self.target_wavelength[i],
                    "nm: ",
                    target_calculated,
                    " != ",
                    self.target_value[i],
                )
                test_pass = False

            if self.target_condition[i] == ">" and target_calculated <= float(
                self.target_value[i]
            ):

                print(
                    "Target not met: at ",
                    self.target_wavelength[i],
                    "nm: ",
                    target_calculated,
                    " <= ",
                    self.target_value[i],
                )
                test_pass = False

            if self.target_condition[i] == "<" and target_calculated > float(
                self.target_value[i]
            ):

                print(
                    "Target not met: at ",
                    self.target_wavelength[i],
                    "nm: ",
                    target_calculated,
                    " >= ",
                    self.target_value[i],
                )
                test_pass = False

        return test_pass

    def perform_optimisation(self, optimisation_type):

        data_optim = self.data_order

        if optimisation_type == "differential_evolution":
            ret = differential_evolution(
                self.merit_function,
                bounds=self.bounds,
            )
        elif optimisation_type == "dual_annealing":
            ret = dual_annealing(
                self.merit_function,
                bounds=self.bounds,
            )
        elif optimisation_type == "minimize":
            ret = minimize(
                self.merit_function,
                x0=self.initial_thicknesses,
                bounds=self.bounds,
                method="Nelder-Mead",
            )

        return ret.x


if __name__ == "__main__":

    # Simple test to see if everything is working
    optim_module = OptimModule("test_optimisation.json")
    print(optim_module.perform_optimisation("minimize"))
