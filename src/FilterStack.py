import ctypes
import numpy as np
import json
import math
import time
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import (
    dual_annealing,
    minimize,
    differential_evolution,
    basinhopping,
    shgo,
    brute,
)


def ignore(msg=None):
    pass


class FilterStack:
    """
    This class represents an optical filter stack and provides methods for
    its optimization.

    The FilterStack class is initialized with an optimization order file and
    optional parameters for message and update queues, as well as logging
    functions. It translates the optimization order file into a format readable
    by the C++ code that performs the actual filter design and optimization and
    interfaces to C++, too.

    The class uses various optimization algorithms from the scipy.optimize
    module, such as dual_annealing, minimize, differential_evolution,
    basinhopping, shgo, and brute.
    """

    def __init__(
        self,
        json_file_path,
        message_queue=None,
        update_queue=None,
        log_func=print,
        log_design_func=ignore,
    ):
        """
        Initializes the FilterStack class.

        This method translates the provided JSON file into a format readable by
        the C++ code, creates the filter stack in C++, and reads in the
        translated JSON file again for easy handling in Python. It also converts
        the target ranges in the filter definition into arrays of individual
        targets. It sets up helper variables for storing states during
        optimization and logging and queue variables needed for the web-based
        GUI.

        Attributes:
        json_file_path (str): The path to the file containing the optimization order.
        message_queue (Queue, optional): A queue for inter-thread communication to send messages to the GUI.
        update_queue (Queue, optional): A queue for inter-thread communication to send updates to the GUI.
        log_func (function, optional): A function for logging general messages. Defaults to the print function.
        log_design_func (function, optional): A function for logging design-specific
        messages. Defaults to an ignore function that does nothing.
        """

        # translate json file to a readable format for the C++ code that does
        # not involve abbreviations that we use for filter design
        json_file_path_cpp = self.translate_order_for_cpp(json_file_path)

        # Create the filter stack in C++
        self.my_filter, self.lib = self.create_filter_in_cpp(json_file_path_cpp)

        # Read in the json file translated for C++ again and restructure a bit
        # for easy handling in python
        with open(json_file_path_cpp) as f:
            self.filter_definition = json.load(f)

        self.filter_definition["structure_thicknesses"] = np.array(
            self.filter_definition["structure_thicknesses"]
        )
        self.bounds = [tuple(x) for x in self.filter_definition["bounds"]]

        (
            self.target_wavelength,
            self.target_polar_angle,
            self.target_azimuthal_angle,
            self.target_weights,
            self.target_value,
            self.target_polarization,
            self.target_condition,
            self.target_tolerance,
            self.target_type,
        ) = self.convert_range_targets(self.filter_definition)

        # Helper variables that store states during optimization
        self.last_log_time = 0
        self.optimum_merit = None
        self.optimum_number = 0
        self.optimum_iteration = None
        self.last_optimum_number = 0
        self.first_zero = True
        self.is_general_core = True

        self.initial_merit = 0
        self.iteration_no = 0
        self.callback_call = 0
        self.previous_layer_positions = np.arange(
            0, np.sum(self.filter_definition["layer_switch_allowed"]), 1
        )

        # Logging and queue variables needed for web-based GUI
        self.log_func = log_func
        self.log_design_func = log_design_func
        self.message_queue = message_queue
        self.update_queue = update_queue

    #####################################################
    ######### Initialization and C++ Interfacing ########
    #####################################################
    def translate_order_for_cpp(self, json_file_path):
        """
        Translates the optimization order for the C++ code.

        This function reads in a JSON file that specifies the optimization order for the filter stack. It translates the structure of the filter stack, including the materials, thicknesses, and optimization parameters, into a format that can be read by the C++ code. If the structure includes periodic assemblies, these are expanded into individual layers. If additional layers are allowed to be added during optimization, these are also included in the translated order. The translated order is then written to a temporary JSON file, which is returned by the function.

        Parameters:
        json_file_path (str): The path to the JSON file that specifies the optimization order.

        Returns:
        str: The path to the temporary JSON file that contains the translated optimization order.
        """

        with open(json_file_path, "r") as json_file_path:
            optimisation_order = json.load(json_file_path)

        updated_optimisation_order = optimisation_order.copy()

        # translate the structure

        structure_materials = optimisation_order["structure_materials"]
        structure_thicknesses = optimisation_order["structure_thicknesses"]
        thickness_opt_allowed = optimisation_order["thickness_opt_allowed"]
        layer_switch_allowed = optimisation_order["layer_switch_allowed"]
        add_layers = optimisation_order["add_layers"]
        bounds = optimisation_order["bounds"]

        updated_structure_materials = []
        updated_structure_thicknesses = []
        updated_thickness_opt_allowed = []
        updated_layer_switch_allowed = []
        updated_bounds = []

        for m_idx, mat in enumerate(structure_materials):

            ## periodic assemblies start with a number

            if mat[0].isdigit():
                parts = mat.split("_")

                nb_layers = int(parts[0])
                layer1_mat = parts[1]
                layer2_mat = parts[2]

                for i in range(nb_layers):
                    updated_structure_materials.append(layer1_mat)
                    updated_structure_materials.append(layer2_mat)
                    updated_structure_thicknesses.append(
                        structure_thicknesses[m_idx][0]
                    )
                    updated_structure_thicknesses.append(
                        structure_thicknesses[m_idx][1]
                    )
                    updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
                    updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
                    updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
                    updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
                    updated_bounds.append(bounds[m_idx])
                    updated_bounds.append(bounds[m_idx])

            else:
                updated_structure_materials.append(mat)
                updated_structure_thicknesses.append(structure_thicknesses[m_idx])
                updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
                updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
                updated_bounds.append(bounds[m_idx])

        updated_optimisation_order["structure_materials"] = updated_structure_materials
        updated_optimisation_order["structure_thicknesses"] = (
            updated_structure_thicknesses
        )
        updated_optimisation_order["thickness_opt_allowed"] = (
            updated_thickness_opt_allowed
        )
        updated_optimisation_order["layer_switch_allowed"] = (
            updated_layer_switch_allowed
        )
        updated_optimisation_order["bounds"] = updated_bounds
        if add_layers:
            for i in range(0, optimisation_order["nb_added_layers"]):
                updated_optimisation_order["structure_materials"].append(
                    optimisation_order["added_materials"][i]
                )
                updated_optimisation_order["structure_thicknesses"].append(
                    0.5
                    * (
                        optimisation_order["added_layer_bounds"][i][0]
                        + optimisation_order["added_layer_bounds"][i][1]
                    )
                )
                updated_optimisation_order["bounds"].append(
                    [
                        optimisation_order["added_layer_bounds"][i][0],
                        optimisation_order["added_layer_bounds"][i][1],
                    ]
                )
                updated_optimisation_order["thickness_opt_allowed"].append(True)
                updated_optimisation_order["layer_switch_allowed"].append(True)

        with open("./temp/temp_cpp_order.json", "w") as f:
            json.dump(updated_optimisation_order, f)

        return "./temp/temp_cpp_order.json"

    def create_filter_in_cpp(self, json_file_path_cpp):
        """
        Creates a filter stack in C++ and returns a pointer to it.

        This function loads a shared library that contains the C++
        implementation of the filter stack. It sets up the argument types and
        return types for the C++ functions that will be called from Python. It
        then calls the createFilterStack function in the C++ library to create a
        new filter stack based on the provided JSON file, and returns a pointer
        to the created filter stack and the loaded library.

        Parameters:
        json_file_path_cpp (str): The path to the JSON file that defines the
        filter stack.

        Returns:
        tuple: A tuple containing a pointer to the created filter stack and the
        loaded library.
        """

        # Link C++ functions and create filter
        lib = ctypes.CDLL("./run_filter_stack.so")

        FilterStack = ctypes.POINTER(ctypes.c_char)

        c_double_array = np.ctypeslib.ndpointer(
            dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"
        )

        c_int_array = np.ctypeslib.ndpointer(
            dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"
        )

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

        lib.change_material_thickness.argtypes = [
            FilterStack,
            c_double_array,
            ctypes.c_size_t,
        ]
        lib.change_material_order.argtypes = [FilterStack, c_int_array, ctypes.c_size_t]
        lib.reset_filter.argtypes = [FilterStack]
        lib.get_material_order.argtypes = [FilterStack]
        lib.get_thicknesses.argtypes = [FilterStack]

        my_filter = lib.createFilterStack(json_file_path_cpp.encode("utf-8"))

        return my_filter, lib

    def convert_range_targets(self, filter_definition):
        """
        Converts the target ranges in the filter definition into arrays of
        individual targets.

        This function takes a filter definition dictionary that includes ranges
        of target wavelengths, polar angles, and azimuthal angles. It converts
        these ranges into arrays of individual targets, each with its own weight
        for the merit function evaluation. If the targets are specified as
        individual values instead of ranges, they are converted into
        single-element arrays.

        Parameters:
        filter_definition (dict): The filter definition, which includes the target ranges.

        Returns:
        tuple: A tuple containing arrays of target wavelengths, polar angles,
        azimuthal angles, weights, values, polarizations, conditions,
        tolerances, and types.
        """

        target_wavelength = np.empty(0)
        target_weights = np.empty(0)
        target_value = np.empty(0)
        target_polarization = np.empty(0)
        target_condition = np.empty(0)
        target_polar_angle = np.empty(0)
        target_azimuthal_angle = np.empty(0)
        target_type = np.empty(0)
        target_tolerance = np.empty(0)

        for target_idx in range(0, np.size(filter_definition["targets_type"])):

            # polar angle treament in case of intervals, each angle is weighted
            # for evaluation with the merit function.

            if isinstance(filter_definition["targets_polar_angle"][target_idx], list):

                polar_angle_entry = np.array(
                    filter_definition["targets_polar_angle"][target_idx]
                )

                interval_polar = np.arange(
                    polar_angle_entry[0],
                    polar_angle_entry[-1] + 1,
                    filter_definition["polar_angle_steps"][target_idx],
                )

                weight_polar = 1 / np.size(interval_polar)

            else:

                interval_polar = [filter_definition["targets_polar_angle"][target_idx]]
                weight_polar = 1

            # azimuthal angle treament in case of intervals, each angle is weighted
            # for evaluation with the merit function.

            if isinstance(
                filter_definition["targets_azimuthal_angle"][target_idx], list
            ):

                azimuthal_entry = np.array(
                    filter_definition["targets_azimuthal_angle"][target_idx]
                )

                interval_azim = np.arange(
                    azimuthal_entry[0],
                    azimuthal_entry[-1] + 1,
                    filter_definition["azimuthal_angle_steps"][target_idx],
                )

                weight_azim = 1 / np.size(interval_azim)

            else:

                interval_azim = [
                    filter_definition["targets_azimuthal_angle"][target_idx]
                ]
                weight_azim = 1

            # wavelength treament in case of intervals, each wavelength is weighted
            # for evaluation with the merit function.

            if isinstance(filter_definition["targets_wavelengths"][target_idx], list):

                wavelength_entry = np.array(
                    filter_definition["targets_wavelengths"][target_idx]
                )

                interval_wvl = np.arange(
                    wavelength_entry[0],
                    wavelength_entry[-1] + 1,
                    filter_definition["wavelength_steps"][target_idx],
                )

                weight_wvl = 1 / np.size(interval_wvl)

            else:

                interval_wvl = [filter_definition["targets_wavelengths"][target_idx]]
                weight_wvl = 1

            for polar_angle in interval_polar:

                # print("polar_angle: ", polar_angle)

                for azim_angle in interval_azim:

                    # print("azim_angle: ", azim_angle)

                    for wvl in interval_wvl:

                        # print("wavelength: ", wvl)

                        target_wavelength = np.append(target_wavelength, wvl)

                        target_polar_angle = np.append(target_polar_angle, polar_angle)

                        target_azimuthal_angle = np.append(
                            target_azimuthal_angle, azim_angle
                        )
                        target_weights = np.append(
                            target_weights, weight_wvl * weight_azim * weight_polar
                        )

                        target_value = np.append(
                            target_value,
                            filter_definition["targets_value"][target_idx],
                        )

                        target_polarization = np.append(
                            target_polarization,
                            filter_definition["targets_polarization"][target_idx],
                        )

                        target_condition = np.append(
                            target_condition,
                            filter_definition["targets_condition"][target_idx],
                        )

                        target_tolerance = np.append(
                            target_tolerance,
                            filter_definition["targets_tolerance"][target_idx],
                        )

                        target_type = np.append(
                            target_type,
                            filter_definition["targets_type"][target_idx],
                        )
        return (
            target_wavelength,
            target_polar_angle,
            target_azimuthal_angle,
            target_weights,
            target_value,
            target_polarization,
            target_condition,
            target_tolerance,
            target_type,
        )

    #####################################################
    ############# Filter Optimization Area ##############
    #####################################################

    def perform_optimisation(
        self, optimisation_type, save_optimized_to_file=True, stop_flag=None
    ):
        """
        Performs the optimization process based on the specified optimization
        type. This function sets up the initial conditions and bounds for the
        optimization, then performs the optimization using the specified method.
        It extracts the optimized thicknesses and layer order from the result
        and applies them to the filter. It then logs the optimization time,
        optimized features, merit value, and number of function evaluations, and
        optionally saves these values to a file.

        Parameters:
        optimisation_type (str): The type of optimization to perform. This can
        be "dual_annealing", "differential_evolution", "basinhopping", "brute",
        "shgo", or "minimize".
        save_optimized_to_file (bool, optional): Whether to save the optimized
        values to a file. Defaults to True.
        stop_flag (function, optional): A function that returns True if the
        optimization should be stopped, False otherwise. Defaults to None.

        Returns:
        numpy array: The optimized features.
        """

        print("running optimisation...")
        self.stop_flag = stop_flag
        start_time = time.time()

        x_initial = []
        bounds = []

        if np.any(self.filter_definition["thickness_opt_allowed"]):
            # Assemble initial values and bounds depending on the layers that
            # were selected for optimization
            x_initial = (
                x_initial
                + self.filter_definition["structure_thicknesses"][
                    self.filter_definition["thickness_opt_allowed"]
                ].tolist()
            )

            # There must be a better way of indexing!
            for i in range(np.shape(self.bounds)[0]):
                if self.filter_definition["thickness_opt_allowed"][i]:
                    bounds.append(self.bounds[i])

        if np.any(self.filter_definition["layer_switch_allowed"]):
            # If layer switching is allowed add the additional parameter that
            # decides the layer positions

            # Evently distribute the initial layer positions over the stack to
            # minimize interference during optimization.

            initial_positions = np.where(
                self.filter_definition["layer_switch_allowed"]
            )[0].tolist()

            x_initial = x_initial + initial_positions

            # Set the first layer positions to the initial positions
            self.previous_layer_positions = initial_positions

            for i in range(np.sum(self.filter_definition["layer_switch_allowed"])):
                bounds.append(
                    (
                        -0.5,
                        np.size(self.filter_definition["structure_thicknesses"]) + 0.5,
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
            with open("./temp/optimized_values.csv", "w") as the_file:
                the_file.write("\n".join(optimized_values))

        return ret.x

    def merit_function(self, features):
        """
        Computes the merit function for the given features. This function takes
        a list of features, extracts the layer thicknesses and positions, and
        computes a merit value based on how well the current configuration of
        the filter stack meets the target conditions. The merit value is a
        measure of the difference between the target and calculated values,
        normalized by the target tolerance and weighted by the target weights.
        The function also checks for stop conditions and raises an exception if
        the optimization should be stopped.

        Parameters:
        features (numpy array): The input array of features.

        Returns:
        float: The computed merit value, normalized by the initial merit value.
        """

        # Declare merit
        merit = 0

        # Extract thickness and layer order from features
        thicknesses, layer_order = self.extract_thickness_and_position_from_features(
            features
        )

        # Loop over targets
        for i in range(0, np.size(self.target_value)):

            # Change material thickness
            if np.any(self.filter_definition["thickness_opt_allowed"]):
                self.lib.change_material_thickness(
                    self.my_filter, thicknesses, int(np.size(thicknesses))
                )

            # Change layer switch allowed
            if np.any(self.filter_definition["layer_switch_allowed"]):
                self.lib.change_material_order(
                    self.my_filter, layer_order, int(np.size(layer_order))
                )

            # Select which core to use
            if self.filter_definition["core_selection"] == "general":
                self.is_general_core = True
            elif self.filter_definition["core_selection"] == "fast":
                self.is_general_core = False
            else:
                if self.target_polarization[i] != "s":
                    self.is_general_core = self.lib.getGeneralMaterialsInStack(
                        [self.my_filter]
                    )
                else:
                    self.is_general_core = False

            # Calculate actual value
            target_calculated = self.lib.calculate_reflection_transmission_absorption(
                self.my_filter,
                self.target_type[i].encode("utf-8"),
                self.target_polarization[i].encode("utf-8"),
                float(self.target_wavelength[i]),
                float(self.target_polar_angle[i]),
                float(self.target_azimuthal_angle[i]),
                self.is_general_core,
            )

            # Calculate merit for equal condition
            if self.target_condition[i] == "=" and target_calculated != float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_weights[i]

            # Calculate merit for larger condition
            if self.target_condition[i] == ">" and target_calculated < float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_weights[i]

            # Calculate merit for smaller condition
            if self.target_condition[i] == "<" and target_calculated > float(
                self.target_value[i]
            ):

                merit += (
                    (target_calculated - self.target_value[i])
                    / float(self.target_tolerance[i])
                ) ** 2 * self.target_weights[i]

        # If merit is zero, the optimization has reached the target
        if merit == 0:

            # Callback once to save the obtained features
            if self.iteration_no == 0:
                raise Exception("Initial features are optimal.")

            if self.first_zero:
                self.first_zero = False
                self.log_func("Optimization has reached merit 0")
                self.callback(features, merit)

            return 0

        else:

            # Set initial merit and normalize to it
            if self.initial_merit == 0:
                self.initial_merit = merit
                self.optimum_merit = merit

            ## callback function
            if self.callback(features, merit):
                raise Exception("Optimization stopped.")

            self.iteration_no += 1

            return merit / self.initial_merit

    def callback(self, x, f):
        """
        Callback function that is called after each optimization step to impose
        further break conditions and log the progress. This function checks if
        the current merit function value `f` is less than the previously best
        found value or equal to zero. If so, it updates the optimum values and
        saves the current structure to a JSON file. It also logs the progress
        every 2 seconds and checks for stop conditions. This is called via the
        merit function and not via scipy's callbacks as this is known to not
        work reliably.

        Parameters:
        x (numpy array): The current solution vector in the optimization process.
        f (float): The current value of the merit function.
        stop_flag (function): A function that returns True if the optimization
        should be stopped, False otherwise.

        Returns:
        bool: True if the optimization should be stopped, False otherwise.
        """

        # Save the current best optimisation values to file
        if f < self.optimum_merit or f == 0:
            temp_json = self.filter_definition.copy()
            if (
                np.any(self.filter_definition["layer_switch_allowed"])
                or self.filter_definition["add_layers"]
            ):
                thicknesses, positions = (
                    self.extract_thickness_and_position_from_features(x)
                )
                current_structure_indices = positions.astype(np.int32)
                current_structure_materials = [
                    self.filter_definition["structure_materials"][el]
                    for el in current_structure_indices
                ]
                current_structure_thicknesses = list(x)
                temp_json["structure_materials"] = current_structure_materials
                temp_json["structure_thicknesses"] = current_structure_thicknesses

                if self.filter_definition["add_layers"]:
                    temp_json["add_layers"] = False
            else:
                temp_json["structure_thicknesses"] = list(x)

            self.optimum_merit = f
            self.optimum_iteration = self.iteration_no
            self.optimum_number += 1

            with open("./temp/current_structure.json", "w") as file:
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

        # If the merit function is close to zero (within the tolerance) or the
        # stop flag is true, stop the optimization.
        if self.stop_flag() or math.isclose(f, 0, abs_tol=1e-6):
            return True
        else:
            return False

    def extract_thickness_and_position_from_features(self, features):
        """
        Extracts layer thicknesses and positions from the given features. This
        function takes a list of features and determines whether each feature
        represents a layer thickness or a layer position based on the
        `thickness_opt_allowed` and `layer_switch_allowed` attributes. It then
        extracts the thicknesses and positions and returns them as separate
        arrays.

        Parameters:
        features (numpy array): The input array of features.

        Returns:
        tuple: A tuple containing two numpy arrays. The first array contains the
        layer thicknesses and the second array contains the layer positions.
        """

        if np.any(self.filter_definition["thickness_opt_allowed"]) and not np.any(
            self.filter_definition["layer_switch_allowed"]
        ):
            # All features are thicknesses
            thicknesses = np.copy(self.filter_definition["structure_thicknesses"])
            thicknesses[np.where(self.filter_definition["thickness_opt_allowed"])] = (
                features
            )
            thicknesses = thicknesses.astype(np.float64)
            layer_order = np.arange(
                0, len(self.filter_definition["structure_thicknesses"])
            )
        elif not np.any(self.filter_definition["thickness_opt_allowed"]) and np.any(
            self.filter_definition["layer_switch_allowed"]
        ):
            # All features are layer positions
            thicknesses = np.copy(self.filter_definition["structure_thicknesses"])

            # Transform the layer positions to a full stack order
            layer_order = self.convert_layer_positions_to_stack_order(features)

            # layer_order = self.allowed_permutations[
            # self.clamp(features[-1], 0, len(self.allowed_permutations) - 1)
            # ].astype(np.int32)
        elif np.any(self.filter_definition["thickness_opt_allowed"]) and np.any(
            self.filter_definition["layer_switch_allowed"]
        ):
            # Some features are thicknesses and some are layer positions
            thicknesses = np.copy(self.filter_definition["structure_thicknesses"])
            thicknesses[np.where(self.filter_definition["thickness_opt_allowed"])] = (
                features[: -1 * np.sum(self.filter_definition["layer_switch_allowed"])]
            )
            thicknesses = thicknesses.astype(np.float64)

            # The hard part is to find the correct way of generating a sensible
            # structure from the feature numbers
            temp_layer_positions = features[
                -1 * np.sum(self.filter_definition["layer_switch_allowed"]) :
            ]

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

    def convert_layer_positions_to_stack_order(self, temp_layer_positions):
        """
        Converts the given layer positions into a full device stack order. This
        function takes an array of layer positions and converts it into a full
        device stack order. It first rounds the layer positions to the nearest
        integer, then clamps them to the available number of layers. It ensures
        that each position is occupied by only one layer. If a layer is allowed
        to move, it removes it from its current position and inserts it at the
        new position.

        Parameters:
        temp_layer_positions (numpy array): The input array of layer positions.

        Returns:
        numpy array: The modified array representing the full device stack order.
        """

        # Now ceil or floor values depening on a threshold value (to deal
        # with the integer character of the features)
        # threshold = 1e-4
        # temp_layer_positions = np.where(temp_layer_positions - self.previous_layer_positions > threshold, np.ceil(temp_layer_positions), np.where(temp_layer_positions - self.previous_layer_positions < threshold, np.floor(temp_layer_positions), temp_layer_positions))
        temp_layer_positions = np.round(temp_layer_positions, 0)

        # Now clamp the integers to the available number of layers
        temp_layer_positions = np.clip(
            temp_layer_positions,
            0,
            len(self.filter_definition["structure_thicknesses"]),
        ).astype(np.int32)

        # If feature 1 has a certain number, the second layer cannot have
        # the same (one position can only be occupied by one layer)
        unique_layer_positions = self.ensure_unique(temp_layer_positions)

        # Set previous layer positions to be the current ones
        self.previous_layer_positions = unique_layer_positions

        # Start with the initial order of layers
        layer_order = list(range(len(self.filter_definition["structure_thicknesses"])))

        # Iterate over each layer
        j = 0
        for i in range(len(self.filter_definition["layer_switch_allowed"])):
            if self.filter_definition["layer_switch_allowed"][i]:
                # If the layer is allowed to move, remove it from its current
                # position and insert it at the new position
                layer_order.remove(i)
                layer_order.insert(unique_layer_positions[j], i)
                j += 1

        return np.array(layer_order).astype(np.int32)

    def ensure_unique(self, arr):
        """
        Ensures all elements in the input array are unique. This function
        iterates over the array. If it finds a duplicate element, it increments
        the value of the element until it is unique.

        Parameters:
        arr (list): The input array.

        Returns:
        list: The modified array with unique elements.
        """

        for i in range(1, len(arr)):
            while arr[i] in arr[:i]:
                arr[i] += 1
        return arr

    #####################################################
    ########## Plotting and calculation area ############
    #####################################################

    def calculate_ar_data(
        self,
        wavelength=None,
        polar_angles=None,
        azimuthal_angles=None,
        target_type=None,
        polarization=None,
        save_figure=False,
        save_data=False,
        web=False,
    ):
        """
        Calculates the angle-resolved (AR)) data for a given set of parameters.

        Parameters:
        wavelength (numpy array, optional): The wavelengths to consider. If not provided, a range is generated based on the filter_definition attribute.
        polar_angles (numpy array, optional): The polar angles to consider. If not provided, a range is generated based on the filter_definition attribute.
        azimuthal_angles (numpy array, optional): The azimuthal angles to consider. Not used in the current implementation.
        target_type (str, optional): The type of target to consider. Not used in the current implementation.
        polarization (str, optional): The type of polarization to consider. Not used in the current implementation.
        save_figure (bool, optional): Whether to save the resulting figure. Defaults to False.
        save_data (bool, optional): Whether to save the resulting data. Defaults to False.
        web (bool, optional): Whether to prepare the data for web usage. Defaults to False.

        Returns:
        The function plots the data to file or display
        """
        if wavelength is None:
            wavelength = np.arange(
                self.filter_definition["wavelengthMin"],
                self.filter_definition["wavelengthMax"] + 1,
                self.filter_definition["wavelengthStep"],
            )
        if polar_angles is None:
            polar_angles = np.arange(
                self.filter_definition["polarAngleMin"],
                self.filter_definition["polarAngleMax"] + 1,
                self.filter_definition["polarAngleStep"],
            )
        if azimuthal_angles is None:
            azimuthal_angles = np.arange(
                self.filter_definition["azimAngleMin"],
                self.filter_definition["azimAngleMax"] + 1,
                self.filter_definition["azimAngleStep"],
            )

        if target_type is None:
            target_type = self.filter_definition["calculation_type"]

        if polarization is None:
            polarization = self.filter_definition["polarization"]

        if self.filter_definition["core_selection"] == "general":
            is_general_core = True
        elif self.filter_definition["core_selection"] == "fast":
            is_general_core = False
        else:
            if polarization:
                is_general_core = self.lib.getGeneralMaterialsInStack([self.my_filter])
            else:
                is_general_core = False

        initial_time = time.time()
        stored_value = pd.DataFrame(columns=polar_angles.astype("U"), index=wavelength)
        for phi in azimuthal_angles:
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

                if web:
                    self.log_func(
                        {
                            "intensity": np.nan_to_num(
                                stored_value.to_numpy(float)
                            ).tolist(),
                            "wavelength": stored_value.index.astype(float).tolist(),
                            "angles": stored_value.columns.astype(float).tolist(),
                            "azimuthal_angle": phi,
                        }
                    )

        self.log_design_func()

        if not web:
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
                plt.savefig("./temp/plot.png", format="png", dpi=300)
                # Save X, Y, Z to csv files
            if save_data:
                header_lines = []
                header_lines.append("Here we can add some meta data to the header \n")

                # Save header lines indicating what the simulation represents
                with open("./temp/value.csv", "w") as the_file:
                    the_file.write("\n".join(header_lines))

                # Save actual data by appending
                stored_value.to_csv("./temp/value.csv", sep=",", header=True, mode="a")

            plt.show()

    #####################################################
    ################### Unit Testing ####################
    #####################################################

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
