# import ctypes
import numpy as np
import json
import itertools
import math
import time

from scipy.optimize import (
    dual_annealing,
    minimize,
    differential_evolution,
    basinhopping,
    shgo,
    brute,
)

from queue import Queue


def ignore(msg=None):
    pass


class OptimModule:
    def __init__(
        self,
        optimisation_order_file,
        my_filter,
        ctypes_lib,
        message_queue=None,
        update_queue=None,
        log_func=print,
        log_design_func=ignore,
    ):

        self.last_log_time = 0
        self.optimum_merit = None
        self.optimum_number = 0
        self.optimum_iteration = None
        self.last_optimum_number = 0
        self.first_zero = True
        self.lib = ctypes_lib
        self.my_filter = my_filter
        self.log_func = log_func
        self.log_design_func = log_design_func
        self.message_queue = message_queue
        self.update_queue = update_queue
        with open(optimisation_order_file) as f:
            self.data_order = json.load(f)
        self.initial_thicknesses = np.array(self.data_order["structure_thicknesses"])
        self.thickness_opt_allowed = np.array(self.data_order["thickness_opt_allowed"])
        self.layer_switch_allowed = np.array(self.data_order["layer_switch_allowed"])
        self.type_entries = np.array(self.data_order["targets_type"])
        self.polarization_entries = np.array(self.data_order["targets_polarization"])
        self.value_entries = np.array(self.data_order["targets_value"])
        self.condition_entries = np.array(self.data_order["targets_condition"])
        self.polar_angle_entries = self.data_order["targets_polar_angle"]
        self.polar_angle_steps = np.array(self.data_order["polar_angle_steps"])
        self.azim_angle_entries = self.data_order["targets_azimuthal_angle"]
        self.azim_angle_steps = np.array(self.data_order["azimuthal_angle_steps"])
        self.wavelength_entries = self.data_order["targets_wavelengths"]
        self.wavelength_step = np.array(self.data_order["wavelength_steps"])
        self.tolerance_entries = np.array(self.data_order["targets_tolerance"])
        self.bounds = [tuple(x) for x in self.data_order["bounds"]]
        self.core_override = None
        if self.data_order["core_override"]:
            self.core_override = self.data_order["core_override"]
            if self.data_order["core_selected"] == "general":
                self.is_general_core = True
            else:
                self.is_general_core = False
        self.target_wavelength = np.empty(0)
        self.target_weights = np.empty(0)
        self.target_value = np.empty(0)
        self.target_polarization = np.empty(0)
        self.target_condition = np.empty(0)
        self.target_polar_angle = np.empty(0)
        self.target_azimuthal_angle = np.empty(0)
        self.target_type = np.empty(0)
        self.target_tolerance = np.empty(0)

        self.initial_merit = 0
        self.iteration_no = 0
        self.callback_call = 0
        self.previous_layer_positions = np.arange(
            0, np.sum(self.layer_switch_allowed), 1
        )

        for target_idx in range(0, np.size(self.type_entries)):

            # polar angle treament in case of intervals, each angle is weighted
            # for evaluation with the merit function.

            if isinstance(self.polar_angle_entries[target_idx], list):

                polar_angle_entry = np.array(self.polar_angle_entries[target_idx])

                interval_polar = np.arange(
                    polar_angle_entry[0],
                    polar_angle_entry[-1] + 1,
                    self.polar_angle_steps[target_idx],
                )

                weight_polar = 1 / np.size(interval_polar)

            else:

                interval_polar = [self.polar_angle_entries[target_idx]]
                weight_polar = 1

            # azimuthal angle treament in case of intervals, each angle is weighted
            # for evaluation with the merit function.

            if isinstance(self.azim_angle_entries[target_idx], list):

                azimuthal_entry = np.array(self.azim_angle_entries[target_idx])

                interval_azim = np.arange(
                    azimuthal_entry[0],
                    azimuthal_entry[-1] + 1,
                    self.azim_angle_steps[target_idx],
                )

                weight_azim = 1 / np.size(interval_azim)

            else:

                interval_azim = [self.azim_angle_entries[target_idx]]
                weight_azim = 1

            # wavelength treament in case of intervals, each wavelength is weighted
            # for evaluation with the merit function.

            if isinstance(self.wavelength_entries[target_idx], list):

                wavelength_entry = np.array(self.wavelength_entries[target_idx])

                interval_wvl = np.arange(
                    wavelength_entry[0],
                    wavelength_entry[-1] + 1,
                    self.wavelength_step[target_idx],
                )

                weight_wvl = 1 / np.size(interval_wvl)

            else:

                interval_wvl = [self.wavelength_entries[target_idx]]
                weight_wvl = 1

            for polar_angle in interval_polar:

                # print("polar_angle: ", polar_angle)

                for azim_angle in interval_azim:

                    # print("azim_angle: ", azim_angle)

                    for wvl in interval_wvl:

                        # print("wavelength: ", wvl)

                        self.target_wavelength = np.append(self.target_wavelength, wvl)

                        self.target_polar_angle = np.append(
                            self.target_polar_angle, polar_angle
                        )

                        self.target_azimuthal_angle = np.append(
                            self.target_azimuthal_angle, azim_angle
                        )
                        self.target_weights = np.append(
                            self.target_weights, weight_wvl * weight_azim * weight_polar
                        )

                        self.target_value = np.append(
                            self.target_value, self.value_entries[target_idx]
                        )

                        self.target_polarization = np.append(
                            self.target_polarization,
                            self.polarization_entries[target_idx],
                        )

                        self.target_condition = np.append(
                            self.target_condition, self.condition_entries[target_idx]
                        )

                        self.target_tolerance = np.append(
                            self.target_tolerance, self.tolerance_entries[target_idx]
                        )

                        self.target_type = np.append(
                            self.target_type, self.type_entries[target_idx]
                        )

        # Check that the targets have been entered correctly
        # print("target_wavelength:", self.target_wavelength)
        # print("target_weights:", self.target_weights)
        # print("target_value:", self.target_value)
        # print("target_polarization:", self.target_polarization)
        # print("target_condition:", self.target_condition)
        # print("target_tolerance:", self.target_tolerance)
        # print("target_type:", self.target_type)
        # print("target_polar_angle:", self.target_polar_angle)
        # print("target_azimuthal_angle:", self.target_azimuthal_angle)

    def ensure_unique(self, arr):
        for i in range(1, len(arr)):
            while arr[i] in arr[:i]:
                arr[i] += 1
        return arr

    def convert_layer_positions_to_stack_order(self, temp_layer_positions):
        """
        Convert for given features that stand for a layer position of a certain
        layer to a full device stack order.
        """

        threshold = 5e-1

        # Now ceil or floor values depening on a threshold value (to deal
        # with the integer character of the features)
        # temp_layer_positions = np.where(temp_layer_positions - self.previous_layer_positions > threshold, np.ceil(temp_layer_positions), np.where(temp_layer_positions - self.previous_layer_positions < threshold, np.floor(temp_layer_positions), temp_layer_positions))
        temp_layer_positions = np.round(temp_layer_positions, 0)

        # Now clamp the integers to the available number of layers
        temp_layer_positions = np.clip(
            temp_layer_positions, 0, len(self.initial_thicknesses)
        ).astype(np.int32)

        # If feature 1 has a certain number, the second layer cannot have
        # the same (one position can only be occupied by one layer)
        unique_layer_positions = self.ensure_unique(temp_layer_positions)

        # Set previous layer positions to be the current ones
        self.previous_layer_positions = unique_layer_positions

        # Start with the initial order of layers
        layer_order = list(range(len(self.initial_thicknesses)))

        # Iterate over each layer
        j = 0
        for i in range(len(self.layer_switch_allowed)):
            if self.layer_switch_allowed[i]:
                # If the layer is allowed to move, remove it from its current
                # position and insert it at the new position
                layer_order.remove(i)
                layer_order.insert(unique_layer_positions[j], i)
                j += 1

        return np.array(layer_order).astype(np.int32)

    def extract_thickness_and_position_from_features(self, features):
        """
        If the layer position is to be optimized the last element of the
        features list is the layer position parameter. If it is not to be
        optimized all features are thicknesses.
        """

        if np.any(self.thickness_opt_allowed) and not np.any(self.layer_switch_allowed):
            # All features are thicknesses
            thicknesses = np.copy(self.initial_thicknesses)
            thicknesses[np.where(self.thickness_opt_allowed)] = features
            thicknesses = thicknesses.astype(np.float64)
            layer_order = np.arange(0, len(self.initial_thicknesses))
        elif not np.any(self.thickness_opt_allowed) and np.any(
            self.layer_switch_allowed
        ):
            # All features are layer positions
            thicknesses = np.copy(self.initial_thicknesses)

            # Transform the layer positions to a full stack order
            layer_order = self.convert_layer_positions_to_stack_order(features)

            # layer_order = self.allowed_permutations[
            # self.clamp(features[-1], 0, len(self.allowed_permutations) - 1)
            # ].astype(np.int32)
        elif np.any(self.thickness_opt_allowed) and np.any(self.layer_switch_allowed):
            # Some features are thicknesses and some are layer positions
            thicknesses = np.copy(self.initial_thicknesses)
            thicknesses[np.where(self.thickness_opt_allowed)] = features[
                : -1 * np.sum(self.layer_switch_allowed)
            ]
            thicknesses = thicknesses.astype(np.float64)

            # The hard part is to find the correct way of generating a sensible
            # structure from the feature numbers
            temp_layer_positions = features[-1 * np.sum(self.layer_switch_allowed) :]

            # Transform the layer positions to a full stack order
            layer_order = self.convert_layer_positions_to_stack_order(
                temp_layer_positions
            )

        else:
            # print(
            # "No optimization has been selected for the thicknesses or layer order."
            # )
            raise ValueError
        return thicknesses.astype(np.float64), layer_order.astype(np.int32)

    def callback(self, x, f, stop_flag):
        """
        Function that is called after each optimization step to impose further
        break conditions.
        """
        # Save the current best optimisation values to a file
        if f < self.optimum_merit or f == 0:

            if np.any(self.layer_switch_allowed) or self.data_order["add_layers"]:
                thicknesses, positions = (
                    self.extract_thickness_and_position_from_features(x)
                )
                current_structure_indices = positions.astype(np.int32)
                current_structure_materials = []
                current_structure_thicknesses = []
                for idx, el in enumerate(current_structure_indices):
                    current_structure_materials.append(
                        self.data_order["structure_materials"][el]
                    )
                    current_structure_thicknesses.append(x[idx])
                temp_json = self.data_order.copy()
                temp_json["structure_materials"] = current_structure_materials
                temp_json["structure_thicknesses"] = current_structure_thicknesses

                if self.data_order["add_layers"]:
                    temp_json["add_layers"] = False

            else:
                current_structure_materials = list(
                    self.data_order["structure_materials"]
                )
                current_structure_thicknesses = list(x)
                temp_json = self.data_order.copy()
                temp_json["structure_thicknesses"] = current_structure_thicknesses

            self.optimum_merit = f
            self.optimum_iteration = self.iteration_no
            self.optimum_number += 1

            with open("current_structure.json", "w") as file:
                json.dump(temp_json, file)

        self.callback_call += 1

        current_time = time.time()
        # Number of seconds to wait for next logging is defined here at 2 seconds
        if current_time - self.last_log_time >= 2:
            self.log_func(
                f"merit | call #{str(self.iteration_no)} : {round(f / self.initial_merit, 9)} {round(f, 2)}"
            )
            self.last_log_time = current_time
            if self.last_optimum_number < self.optimum_number:
                self.log_func(
                    f"New optimum on call #{self.optimum_iteration} : {round(self.optimum_merit / self.initial_merit, 9)} {round(self.optimum_merit, 2)}"
                )
                self.log_design_func()
                self.last_optimum_number = self.optimum_number

        # If the merit function is close to zero (within the tolerance), stop
        # the optimization (It probably makes sense here to assume a lower
        # tolerance than for the minimize but this has to be adjusted
        # empirically).
        # If the stop flag is true, the optimization is stopped.

        if self.stop_flag() or math.isclose(f, 0, abs_tol=1e-6):
            return True
        else:
            return False

    def merit_function(self, features):
        """
        Function that adds up all the computed values and computes a merit from the
        results
        """
        merit = 0
        thicknesses, layer_order = self.extract_thickness_and_position_from_features(
            features
        )
        # print(layer_order)
        # print(layer_order)

        # Assemble the whole thicknesses and layer positions from the delivered
        # features. They are set up as follows:
        # [thickness_of_layer_to_optimize_1, .._2, .._3,
        # layer_position_argument]

        for i in range(0, np.size(self.target_value)):

            if np.any(self.thickness_opt_allowed):
                self.lib.change_material_thickness(
                    self.my_filter, thicknesses, int(np.size(thicknesses))
                )
            if np.any(self.layer_switch_allowed):
                self.lib.change_material_order(
                    self.my_filter, layer_order, int(np.size(layer_order))
                )

            if self.core_override == None:
                if self.target_polarization[i] != "s":
                    self.is_general_core = self.lib.getGeneralMaterialsInStack(
                        [self.my_filter]
                    )
                else:
                    self.is_general_core = False
            elif self.core_override != None:
                if self.core_override == "general":
                    self.is_general_core = True
                if self.core_override == "fast":
                    self.is_general_core = False

            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                float(self.target_wavelength[i]),
                float(self.target_polar_angle[i]),
                float(self.target_azimuthal_angle[i]),
                self.is_general_core,
            )

            # print(self.target_wavelength[i], target_calculated)

            if self.target_condition[i] == "=" and target_calculated != float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_weights[i]

            if self.target_condition[i] == ">" and target_calculated < float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_weights[i]

            if self.target_condition[i] == "<" and target_calculated > float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_weights[i]

        if merit != 0:

            # Set initial merit and normalize to it
            if self.initial_merit == 0:
                self.initial_merit = merit
                self.optimum_merit = merit

            ## callback function
            if self.callback(features, merit, False):
                raise Exception("Optimization stopped.")

            self.iteration_no += 1

            return merit / self.initial_merit

        else:

            # Callback once to save the obtained features
            if self.iteration_no == 0:
                raise Exception("Initial features are optimal.")

            if self.first_zero:
                self.first_zero = False
                self.log_func("Optimization has reached merit 0")
                self.callback(features, merit, False)

            return 0

    # check targets is used in unit testing
    def check_targets(self, thicknesses):
        """
        Function that returns true if all targets are met, for unit testing only
        """
        test_pass = True

        for i in range(0, np.size(self.target_value)):

            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                self.target_wavelength[i],
                self.target_polar_angle[i],
                self.target_azimuthal_angle[i],
                self.is_general_core,
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

    def perform_optimisation(
        self, optimisation_type, save_optimized_to_file=True, stop_flag=None
    ):
        self.stop_flag = stop_flag
        start_time = time.time()

        x_initial = []
        bounds = []

        if np.any(self.thickness_opt_allowed):
            # Assemble initial values and bounds depending on the layers that
            # were selected for optimization
            x_initial = (
                x_initial
                + self.initial_thicknesses[self.thickness_opt_allowed].tolist()
            )

            # There must be a better way of indexing!
            for i in range(np.shape(self.bounds)[0]):
                if self.thickness_opt_allowed[i]:
                    bounds.append(self.bounds[i])

        if np.any(self.layer_switch_allowed):
            # If layer switching is allowed add the additional parameter that
            # decides the layer positions

            # Evently distribute the initial layer positions over the stack to
            # minimize interference during optimization.

            initial_positions = np.where(self.layer_switch_allowed)[0].tolist()

            x_initial = x_initial + initial_positions

            # Set the first layer positions to the initial positions
            self.previous_layer_positions = initial_positions

            for i in range(np.sum(self.layer_switch_allowed)):
                bounds.append(
                    (
                        -0.5,
                        np.size(self.initial_thicknesses) + 0.5,
                    )
                )

        ret = 0

        # With scipy we cannot do integer optimization
        if optimisation_type == "dual_annealing":
            ret = dual_annealing(
                self.merit_function,
                bounds=bounds,
            )
        elif optimisation_type == "differential_evolution":
            ret = differential_evolution(
                self.merit_function,
                bounds=bounds,
                # callback = self.callback_func_advanced
            )
        elif optimisation_type == "basinhopping":
            # This algorithm does
            # 1. a random perturbation of the features
            # 2. then a scipy.minimize
            # 3. accepts of rejects the new optimum value
            ret = basinhopping(
                self.merit_function,
                x0=x_initial,
                minimizer_kwargs={
                    "method": "Nelder-Mead",
                },
            )

        elif optimisation_type == "brute":
            # Brute force optimization: only sensible for a low number of
            # features (e.g., 3). Depending on Ns, the number of points to
            # try is established  (Ns^(#features) = number of iterations)
            ret = brute(
                self.merit_function,
                ranges=bounds,
                Ns=2,
                # maxiter = 50000,
            )
        elif optimisation_type == "shgo":
            # Doesn't really start
            ret = shgo(
                self.merit_function,
                bounds=bounds,
                # maxiter = 50000,
            )
        elif optimisation_type == "minimize":

            ret = minimize(
                self.merit_function,
                x0=x_initial,
                bounds=bounds,
                method="Nelder-Mead",
                # the below values for xatol and fatol were found to prevent the function
                # from overoptimising
                options={"xatol": 1e-1, "fatol": 1e-1},
            )

        thicknesses, layer_order = self.extract_thickness_and_position_from_features(
            ret.x
        )
        self.lib.change_material_thickness(
            self.my_filter, thicknesses, int(np.size(thicknesses))
        )
        self.lib.change_material_order(
            self.my_filter, layer_order, int(np.size(layer_order))
        )

        # Set initial merit back to zero
        self.initial_merit = 0
        self.iteration_no = 0
        self.callback_call = 0

        self.log_func("Optimization time: ", time.time() - start_time, "s")
        self.log_func("Optimized features: ", ret.x)
        self.log_func("Optimized merit value: ", ret.fun)
        self.log_func("Number of function evaluations: ", ret.nfev)

        if save_optimized_to_file:
            optimized_values = []
            optimized_values.append(
                "Optimization time: " + str(time.time() - start_time) + " s \n"
            )
            optimized_values.append(
                "Number of function evaluations: " + str(ret.nfev) + "\n"
            )
            optimized_values.append("Optimized merit value: " + str(ret.fun) + "\n")
            optimized_values.append("Optimized features: " + str(ret.x) + "\n")
            with open("optimized_values.csv", "w") as the_file:
                the_file.write("\n".join(optimized_values))

        return ret.x
