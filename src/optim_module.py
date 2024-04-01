# import ctypes
import numpy as np
import json
import itertools
import math
import time

from scipy.optimize import dual_annealing, minimize, differential_evolution, basinhopping, shgo


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
        self.type_entries = np.array(self.data_order["targets_type"])
        self.polarization_entries = np.array(self.data_order["targets_polarization"])
        self.value_entries = np.array(self.data_order["targets_value"])
        self.condition_entries = np.array(self.data_order["targets_condition"])
        self.polar_angle_entries = np.array(self.data_order["targets_polar_angle"])
        self.azim_angle_entries = np.array(self.data_order["targets_azimuthal_angle"])
        self.wavelength_entries = self.data_order["targets_wavelengths"]
        self.wavelength_step = np.array(self.data_order["wavelength_steps"])
        self.tolerance_entries = np.array(self.data_order["targets_tolerance"])
        self.bounds = [tuple(x) for x in self.data_order["bounds"]]

        self.target_wavelength = np.empty(0)
        self.target_wavelength_weights = np.empty(0)
        self.target_value = np.empty(0)
        self.target_polarization = np.empty(0)
        self.target_condition = np.empty(0)
        self.target_polar_angle = np.empty(0)
        self.target_azimuthal_angle = np.empty(0)
        self.target_type = np.empty(0)
        self.target_tolerance = np.empty(0)

        self.initial_merit = 0

        # wavelength treament in case of intervals, each wavelength is weighted
        # for evaluation with the merit function.

        for w_idx, wavelength_entry_unchecked in enumerate(self.wavelength_entries):

            if isinstance(wavelength_entry_unchecked, list):

                wavelength_entry = np.array(wavelength_entry_unchecked)

                interval = np.arange(
                    wavelength_entry[0],
                    wavelength_entry[-1] + 1,
                    self.wavelength_step[w_idx],
                )

                weight = 1 / np.size(interval)

                for wavelength in interval:

                    self.target_wavelength = np.append(
                        self.target_wavelength, wavelength
                    )

                    self.target_wavelength_weights = np.append(
                        self.target_wavelength_weights, weight
                    )

                    self.target_value = np.append(
                        self.target_value, self.value_entries[w_idx]
                    )

                    self.target_polarization = np.append(
                        self.target_polarization, self.polarization_entries[w_idx]
                    )

                    self.target_condition = np.append(
                        self.target_condition, self.condition_entries[w_idx]
                    )

                    self.target_tolerance = np.append(
                        self.target_tolerance, self.tolerance_entries[w_idx]
                    )

                    self.target_type = np.append(
                        self.target_type, self.type_entries[w_idx]
                    )

                    self.target_polar_angle = np.append(
                        self.target_polar_angle, self.polar_angle_entries[w_idx]
                    )

                    self.target_azimuthal_angle = np.append(
                        self.target_azimuthal_angle, self.azim_angle_entries[w_idx]
                    )

            else:

                self.target_wavelength = np.append(
                    self.target_wavelength, wavelength_entry_unchecked
                )
                self.target_wavelength_weights = np.append(
                    self.target_wavelength_weights, 1
                )
                self.target_value = np.append(
                    self.target_value, self.value_entries[w_idx]
                )
                self.target_polarization = np.append(
                    self.target_polarization, self.polarization_entries[w_idx]
                )
                self.target_condition = np.append(
                    self.target_condition, self.condition_entries[w_idx]
                )
                self.target_tolerance = np.append(
                    self.target_tolerance, self.tolerance_entries[w_idx]
                )
                self.target_type = np.append(self.target_type, self.type_entries[w_idx])
                self.target_polar_angle = np.append(
                    self.target_polar_angle, self.polar_angle_entries[w_idx]
                )
                self.target_azimuthal_angle = np.append(
                    self.target_azimuthal_angle, self.azim_angle_entries[w_idx]
                )

        if np.any(self.layer_switch_allowed):
            self.allowed_permutations = self.compute_permutations(
                np.arange(0, len(self.initial_thicknesses)),
                np.where(self.layer_switch_allowed)[0],
            )

    def compute_permutations(self, arr, movable_indices):
        """
        Compute all available permutations for a given stack, given the movable
        layer indices. This is used to permute layers based on a single integer
        value (single feature that needs to be optimized).
        """
        # Remove the movable elements from the array
        movable_elements = arr[movable_indices]
        new_arr = np.delete(arr, movable_indices)

        # Generate all possible positions for the movable elements
        positions = list(itertools.permutations(range(len(arr)), len(movable_elements)))

        # Insert the movable elements into all possible positions
        result = []
        for pos in positions:
            pos = np.array(pos)
            temp_arr = new_arr.copy()
            for idx, element in zip(
                pos[np.argsort(pos)], movable_elements[np.argsort(pos)]
            ):
                temp_arr = np.insert(temp_arr, idx, element)
            result.append(temp_arr)

        return np.array(result)

    def clamp(self, n, minn, maxn):
        """
        Clamp integer values to the maximum number of available positions in the
        stack (given by the number of layers to the power of available movable
        layers).
        """
        return max(min(maxn, int(n)), minn)

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
            thicknesses = np.copy(self.initial_thicknesses)

            layer_order = self.allowed_permutations[
                self.clamp(features[-1], 0, len(self.allowed_permutations) - 1)
            ].astype(np.int32)
        elif np.any(self.thickness_opt_allowed) and np.any(self.layer_switch_allowed):
            thicknesses = np.copy(self.initial_thicknesses)
            thicknesses[np.where(self.thickness_opt_allowed)] = features[:-1]
            thicknesses = thicknesses.astype(np.float64)

            layer_order = self.allowed_permutations[
                self.clamp(features[-1], 0, len(self.allowed_permutations) - 1)
            ].astype(np.int32)
        else:
            print(
                "No optimization has been selected for the thicknesses or layer order."
            )
            raise ValueError
        return thicknesses.astype(np.float64), layer_order.astype(np.int32)

    def callback_func_advanced(self, x, f, context):
        """
        Function that is called after each optimization step to impose further
        break conditions for the non trivial optimization methods.
        """
        # If the merit function is close to zero (within the tolerance), stop
        # the optimization (It probably makes sense here to assume a lower
        # tolerance than for the minimize but this has to be adjusted
        # empirically).
        if math.isclose(f, 0, abs_tol = 1e-5):
            return True
        # print("callback: ", x, f, context)

    def callback_func_minimize(self, intermediate_result):
        """
        Function that is called after each optimization step to impose further
        break conditions for scipy.minimize.
        """
        # If the merit function is close to zero (within the tolerance), stop
        # the optimization.
        if math.isclose(intermediate_result.fun, 0, abs_tol = 1e-9):
            print("Optimization finished because merit reached threshold.")
            raise StopIteration
        # print("callback: ", intermediate_result.fun)

    def merit_function(self, features):
        """
        Function that adds up all the computed values and computes a merit from the
        results
        """
        merit = 0
        # print("Layer order no.: " + str(features))

        thicknesses, layer_order = self.extract_thickness_and_position_from_features(
            features
        )

        # Assemble the whole thicknesses and layer positions from the delivered
        # features. They are set up as follows:
        # [thickness_of_layer_to_optimize_1, .._2, .._3,
        # layer_position_argument]

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
                    self.my_filter, thicknesses, int(np.size(thicknesses))
                )
            if np.any(self.layer_switch_allowed):
                self.lib.change_material_order(
                    self.my_filter, layer_order, int(np.size(layer_order))
                )

            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                float(self.target_wavelength[i]),
                float(self.target_polar_angle[i]),
                float(self.target_azimuthal_angle[i]),
            )

            # print(self.target_wavelength[i], target_calculated)

            if self.target_condition[i] == "=" and target_calculated != float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_wavelength_weights[i]

            if self.target_condition[i] == ">" and target_calculated < float(
                self.target_value[i]
            ):
                # print(
                #     self.target_wavelength[i],
                #     "> : ",
                #     target_calculated,
                #     self.target_value[i],
                # )
                # print(
                #     self.target_wavelength[i],
                #     "> : ",
                #     target_calculated,
                #     float(self.target_value[i]),
                # )

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_wavelength_weights[i]

            if self.target_condition[i] == "<" and target_calculated > float(
                self.target_value[i]
            ):
                # print(
                #     self.target_wavelength[i],
                #     "< : ",
                #     target_calculated,
                #     float(self.target_value[i]),
                # )
                # print(
                #     self.target_wavelength[i],
                #     "< : ",
                #     target_calculated,
                #     float(self.target_value[i]),
                # )

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_wavelength_weights[i]


        # Set initial merit and normalize to it
        if self.initial_merit == 0:
            self.initial_merit = merit

        print("merit: ", merit/self.initial_merit)
        # print("thicknesses: ", thicknesses)
        # print("layer_order: ", layer_order)
        # print("merit: ", merit)
        # print("thicknesses: ", features)

        return merit/self.initial_merit

    # check targets is used in unit testing
    def check_targets(self, thicknesses):
        """
        Function that returns true if all targets are met, for unit testing only
        """
        test_pass = True

        for i in range(0, np.size(self.target_value)):

            # print("target_wavelength:", self.target_wavelength[i])
            # print("target_wavelength_weights:", self.target_wavelength_weights[i])
            # print("target_value:", self.target_value[i])
            # print("target_polarization:", self.target_polarization[i])
            # print("target_condition:", self.target_condition[i])
            # print("target_tolerance:", self.target_tolerance[i])
            # print("target_type:", self.target_type[i])
            # print("target_polar_angle:", self.target_polar_angle[i])
            # print("target_azimuthal_angle:", self.target_azimuthal_angle[i])

            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                self.target_wavelength[i],
                self.target_polar_angle[i],
                self.target_azimuthal_angle[i],
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
            x_initial.append(1)
            bounds.append(
                (
                    -0.5,
                    np.size(self.initial_thicknesses)
                    ** np.sum(self.layer_switch_allowed)
                    + 0.5,
                )
            )

        ret = 0

        # With scipy we cannot do integer optimization
        if optimisation_type == "differential_evolution":
            ret = differential_evolution(
                self.merit_function,
                bounds=bounds,
                callback = self.callback_func_advanced
            )
        elif optimisation_type == "dual_annealing":
            ret = dual_annealing(
                self.merit_function,
                bounds=bounds,
                callback = self.callback_func_advanced,
            )
        elif optimisation_type == "basinhopping":
            ret = basinhopping(
                self.merit_function,
                bounds=bounds,
                callback = self.callback_func_advanced,
            )
        elif optimisation_type == "direct":
            ret = shgo(
                self.merit_function,
                bounds=bounds,
                callback = self.callback_func_advanced,
            )
        elif optimisation_type == "minimize":

            ret = minimize(
                self.merit_function,
                x0=x_initial,
                bounds=bounds,
                method="Nelder-Mead",
                # tol=1e-1,
                callback = self.callback_func_minimize,
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

        print("Optimization time: ", time.time() - start_time, "s")
        print("Optimized features: ", ret.x)

        return ret.x
