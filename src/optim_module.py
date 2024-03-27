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

        self.initial_thicknesses = np.array(self.data_order["structure_thicknesses"])[
            ::-1
        ]
        self.thickness_opt_allowed = np.array(self.data_order["thickness_opt_allowed"])
        self.layer_switch_allowed = np.array(self.data_order["layer_switch_allowed"])
        self.target_type = np.array(self.data_order["targets_type"])
        self.target_polarization = np.array(self.data_order["targets_polarization"])
        self.target_value = np.array(self.data_order["targets_value"])
        self.target_condition = np.array(self.data_order["targets_condition"])
        self.target_angle = np.array(self.data_order["targets_angle"])
        self.target_wavelength = np.array(self.data_order["targets_wavelengths"])
        self.target_tolerance = np.array(self.data_order["targets_tolerance"])
        self.bounds = [tuple(x) for x in self.data_order["bounds"]]

    def merit_function(self, features):
        """
        Function that adds up all the computed values and computes a merit from the
        results
        """
        merit = 0

        # Assemble the whole thicknesses and layer positions from the delivered
        # features. They are set up as follows:
        # [thickness_of_layer_to_optimize_1, .._2, .._3,
        # layer_position_argument]
        if np.any(self.thickness_opt_allowed) and not np.any(self.layer_switch_allowed):
            # All features are thicknesses
            thicknesses = np.copy(self.initial_thicknesses)
            thicknesses[np.where(self.thickness_opt_allowed)] = features
            thicknesses = thicknesses.astype(np.float64)
        else:
            print("This feature was not yet implemented")
            pass

        for i in range(0, np.size(self.target_value)):

            # thicknesses = features[: len(features) / 2].astype(np.float64)
            # layer_positions = features[len(features) / 2 :].astype(np.int32)
            #
            # self.lib.change_material_order(
            # self.my_filter, layer_positions, int(np.size(features) / 2)
            # )
            # self.lib.change_material_thickness(
            # self.my_filter, thicknesses, int(np.size(features) / 2)
            # )
            if np.any(self.thickness_opt_allowed):
                self.lib.change_material_thickness(
                    self.my_filter, thicknesses, np.size(thicknesses)
                )
            if np.any(self.layer_switch_allowed):
                print("This feature has not yet been implemented")
                # self.lib.change_material_order(
                # self.my_filter, layer_positions, int(np.size(features) / 2)
                # )

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
        print("thicknesses: ", features)
        return merit

    # check targets is used in unit testing
    def check_targets(self, thicknesses):
        """
        Function that returns true if all targets are met, for unit testing only
        """
        test_pass = True

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

        if not np.any(self.layer_switch_allowed):
            # Assemble initial values and bounds depending on the layers that
            # were selected for optimization
            x_initial = self.initial_thicknesses[self.thickness_opt_allowed]

            # There must be a better way of indexing!
            bounds = []
            for i in range(np.shape(self.bounds)[0]):
                if self.thickness_opt_allowed[i]:
                    bounds.append(self.bounds[i])
            bounds = np.array(bounds)

        if optimisation_type == "differential_evolution":
            ret = differential_evolution(
                self.merit_function,
                bounds=bounds,
            )
        elif optimisation_type == "dual_annealing":
            ret = dual_annealing(
                self.merit_function,
                bounds=bounds,
            )
        elif optimisation_type == "minimize":
            ret = minimize(
                self.merit_function,
                x0=x_initial,
                bounds=bounds,
                method="Nelder-Mead",
            )

        return ret.x
