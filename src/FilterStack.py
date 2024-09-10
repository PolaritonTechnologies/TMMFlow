# import eventlet

# eventlet.monkey_patch()
import ctypes
import json
import os
import time
import copy
from datetime import datetime
from functools import partial
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

from redis import Redis

from scipy.optimize import (
    #  dual_annealing,
    minimize,
    differential_evolution,
    basinhopping,
    shgo,
    brute,
)

from optimization import (
    dual_annealing,
    gradient_descent,
    gradient,
    particle_swarm,
    grid_search,
)


def ignore(msg=None):
    pass


def log(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message_with_timestamp = f"{timestamp}: {message}"
    print(message_with_timestamp)
    with open("log.txt", "a") as log_file:
        log_file.write(f"{message_with_timestamp}\n")


def convert_numpy_to_list(data):
    if isinstance(data, dict):
        return {key: convert_numpy_to_list(value) for key, value in data.items()}
    elif isinstance(data, np.ndarray):
        return data.tolist()
    elif isinstance(data, list):
        return [convert_numpy_to_list(item) for item in data]
    else:
        return data


def optimization_function(
    optimization_method, latest_design, redis_key, additional_opt_parameters=[]
):
    my_filter = FilterStack(latest_design)
    # database_entry = select_latest_optimization(job_id)
    ret = my_filter.perform_optimisation(
        optimization_method, redis_key, additional_opt_parameters
    )
    return ret
    # db.session.commit()  # Commit changes at the end of the operation


class FilterStack:

    def __init__(
        self,
        my_filter_dict=None,
        my_filter_path=None,
        current_structure="current_structure",
        log_func=log,
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
        json_or_file_path (str): The path to the JSON file that defines the filter stack or the json stack.
        message_queue (Queue, optional): A queue for inter-thread communication to send messages to the GUI.
        update_queue (Queue, optional): A queue for inter-thread communication to send updates to the GUI.
        log_func (function, optional): A function for logging general messages. Defaults to the print function.
        log_design_func (function, optional): A function for logging design-specific
        messages. Defaults to an ignore function that does nothing.
        """
        if my_filter_dict is not None:
            filter_definition_by_user = json.loads(json.dumps(my_filter_dict))
        else:
            with open(my_filter_path, "r") as json_file:
                filter_definition_by_user = json.load(json_file)

        # Initial structure
        self.current_structure = current_structure
        updated_cpp_order = self.translate_order_for_cpp(filter_definition_by_user)

        # Read in the translated C++ again and restructure a bit
        # for easy handling in python
        self.filter_definition = updated_cpp_order
        self.initial_filter_definition = copy.copy(self.filter_definition)

        self.initial_structure_thicknesses = np.copy(
            np.array(self.filter_definition["structure_thicknesses"])
        )
        self.initial_incoherent_order = np.copy(
            np.array(self.filter_definition["incoherent"])
        )

        self.initial_structure_materials = np.copy(
            np.array(self.filter_definition["structure_materials"])
        )

        self.filter_definition["structure_thicknesses"] = np.array(
            self.filter_definition["structure_thicknesses"]
        )
        self.bounds = [tuple(x) for x in self.filter_definition["bounds"]]

        self.layer_order = np.arange(
            0, len(self.filter_definition["structure_thicknesses"]), 1, dtype=np.int32
        )

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
            self.target_arithmetic,
        ) = self.convert_range_targets(self.filter_definition)

        # Helper variables that store states during optimization
        self.last_log_time = 0
        self.optimum_merit = None
        self.optimum_x = None
        self.optimum_number = 0
        self.optimum_iteration = None
        self.last_optimum_number = 0
        self.first_zero = True
        # self.is_general_core = True
        self.stop_flag = False

        self.initial_merit = 0
        self.last_merit = 0
        self.iteration_no = 0
        self.callback_call = 0
        self.previous_layer_positions = np.arange(
            0, np.sum(self.filter_definition["layer_switch_allowed"]), 1
        )

        # Store last calculated data in a pandas dataframe (defaults to None)
        self.stored_data = None

        # Logging and queue variables needed for web-based GUI
        self.log_func = log_func
        self.log_design_func = log_design_func

    #####################################################
    ######### Initialization and C++ Interfacing ########
    #####################################################
    def translate_order_for_cpp(self, optimisation_order):
        """
        Translates the optimization order for the C++ code.

        This function reads in a JSON file that specifies the optimization order
        for the filter stack. It translates the structure of the filter stack,
        including the materials, thicknesses, and optimization parameters, into
        a format that can be read by the C++ code. If the structure includes
        periodic assemblies, these are expanded into individual layers. If
        additional layers are allowed to be added during optimization, these are
        also included in the translated order. The translated order is then
        written to a temporary JSON file, which is returned by the function.

        Parameters:
        json_file_path (str): The path to the JSON file that specifies the optimization order.

        Returns:
        str: The path to the temporary JSON file that contains the translated optimization order.
        """

        updated_optimisation_order = optimisation_order.copy()

        # translate the structure

        structure_materials = optimisation_order["structure_materials"]
        incoherent = optimisation_order["incoherent"]
        structure_thicknesses = optimisation_order["structure_thicknesses"]
        thickness_opt_allowed = optimisation_order["thickness_opt_allowed"]
        layer_switch_allowed = optimisation_order["layer_switch_allowed"]
        # add_layers = optimisation_order["add_layers"]
        bounds = optimisation_order["bounds"]

        updated_structure_materials = []
        updated_incoherent = []
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
                    updated_incoherent.append(incoherent[m_idx])
                    updated_incoherent.append(incoherent[m_idx])
                    updated_bounds.append(bounds[m_idx])
                    updated_bounds.append(bounds[m_idx])

            else:
                updated_structure_materials.append(mat)
                updated_structure_thicknesses.append(structure_thicknesses[m_idx])
                updated_thickness_opt_allowed.append(thickness_opt_allowed[m_idx])
                updated_layer_switch_allowed.append(layer_switch_allowed[m_idx])
                updated_incoherent.append(incoherent[m_idx])
                updated_bounds.append(bounds[m_idx])

        updated_optimisation_order["structure_materials"] = updated_structure_materials
        updated_optimisation_order["incoherent"] = updated_incoherent
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
        """
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
        """

        return updated_optimisation_order

    def create_filter_in_cpp(self):
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
        lib = ctypes.CDLL(
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)), "./cpp/interface.so"
            )
        )

        c_string_array = ctypes.POINTER(ctypes.c_char_p)

        c_double_array = np.ctypeslib.ndpointer(
            dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"
        )

        c_int_array = np.ctypeslib.ndpointer(
            dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"
        )

        lib.createFilterStack.argtypes = [ctypes.c_char_p]
        lib.createFilterStack.restype = c_string_array

        lib.destroyFilterStack.argtypes = [c_string_array]

        lib.calculate_reflection_transmission_absorption_para.argtypes = [
            c_string_array,
            ctypes.c_char_p,
            ctypes.c_double,  # polarization
            c_double_array,
            ctypes.c_size_t,
            c_double_array,
            ctypes.c_size_t,
            c_double_array,
            ctypes.c_size_t,
        ]
        lib.calculate_reflection_transmission_absorption_para.restype = ctypes.c_char_p

        lib.change_material_thickness.argtypes = [
            c_string_array,
            c_double_array,
            ctypes.c_size_t,
        ]
        lib.change_material_order.argtypes = [
            c_string_array,
            c_int_array,
            ctypes.c_size_t,
        ]

        lib.calculate_merit.argtypes = [
            c_string_array,
            c_string_array,
            c_double_array,
            c_double_array,
            c_double_array,
            c_double_array,
            c_double_array,
            c_double_array,
            c_string_array,
            c_double_array,
            c_string_array,
            ctypes.c_size_t,
        ]
        lib.calculate_merit.restype = ctypes.c_double

        lib.initialise_optimization.argtypes = [
            c_string_array,
            c_double_array,
            ctypes.c_size_t,
        ]

        lib.reset_filter.argtypes = [c_string_array]
        lib.get_material_order.argtypes = [c_string_array]
        lib.get_thicknesses.argtypes = [c_string_array]

        # This here is not the most elegant way but currently the only way of
        # creating a new C++ filter from an existing python filter: First create
        # it with the previous thicknesses and layer orders and then reorder and
        # modify thicknesses.
        translated_order = self.translate_order_for_cpp(self.initial_filter_definition)
        my_filter = lib.createFilterStack(json.dumps(translated_order).encode("utf-8"))

        lib.change_material_order(
            my_filter, self.layer_order, int(np.size(self.layer_order))
        )
        lib.change_material_thickness(
            my_filter,
            self.filter_definition["structure_thicknesses"],
            int(np.size(self.filter_definition["structure_thicknesses"])),
        )

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
        target_arithmetic = np.empty(0, dtype="U20")

        old_index = []
        num_elements = 1
        index = 1

        if not np.size(filter_definition["targets_type"]) == 0:
            for target_idx in range(0, np.size(filter_definition["targets_type"])):

                # polar angle treament in case of intervals, each angle is weighted
                # for evaluation with the merit function.

                if isinstance(
                    filter_definition["targets_polar_angle"][target_idx], list
                ):

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

                    interval_polar = [
                        filter_definition["targets_polar_angle"][target_idx]
                    ]
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

                if isinstance(
                    filter_definition["targets_wavelengths"][target_idx], list
                ):

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

                    interval_wvl = [
                        filter_definition["targets_wavelengths"][target_idx]
                    ]
                    weight_wvl = 1

                old_index.append(num_elements)
                num_elements += (
                    len(interval_polar) * len(interval_azim) * len(interval_wvl)
                )

                for polar_angle in interval_polar:

                    for azim_angle in interval_azim:

                        for wvl in interval_wvl:

                            target_wavelength = np.append(target_wavelength, wvl)

                            target_polar_angle = np.append(
                                target_polar_angle, polar_angle
                            )

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

                            if (
                                filter_definition["targets_arithmetic"][target_idx]
                                == ""
                            ):
                                target_arithmetic = np.append(
                                    target_arithmetic,
                                    filter_definition["targets_arithmetic"][target_idx],
                                )
                            else:
                                target_arithmetic = np.append(target_arithmetic, index)

                            index += 1

            # Now actually fill the target_arithmetic with the contents from the initial array to allow for arithmetics
            unexpanded_targets_arithmetic = np.empty(
                np.size(filter_definition["targets_arithmetic"]), dtype="U20"
            )

            # Iterate over the list and replace the numbers
            for i in range(len(unexpanded_targets_arithmetic)):
                # Split the string into components
                components = re.split(
                    "([-+/*])", filter_definition["targets_arithmetic"][i]
                )
                # Replace the numbers with the corresponding values from old_index
                new_components = [
                    str(old_index[int(comp) - 1]) if comp.isdigit() else comp
                    for comp in components
                ]
                # Join the components back together
                unexpanded_targets_arithmetic[i] = "".join(new_components)

            target_arithmetic[np.array(old_index) - 1] = unexpanded_targets_arithmetic
        else:
            log("No optimization targets were defined in the .json.")

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
            target_arithmetic,
        )

    def reset(self):
        """
        Reset filter to initial thickness and order
        """
        self.filter_definition["structure_thicknesses"] = (
            self.initial_structure_thicknesses
        )

        # self.lib.change_material_thickness(
        #     self.my_filter,
        #     self.initial_structure_thicknesses,
        #     int(np.size(self.initial_structure_thicknesses)),
        # )

        self.layer_order = np.arange(
            0, int(np.size(self.initial_structure_thicknesses)), 1, dtype=np.int32
        )

        self.filter_definition["structure_materials"] = [
            self.initial_structure_materials[el] for el in self.layer_order
        ]

        # self.lib.change_material_order(
        #     self.my_filter, self.layer_order, int(np.size(self.layer_order))
        # )

    def return_current_design_as_json(self):
        temp_json = self.filter_definition.copy()

        # Convert to a python list
        temp_json["structure_thicknesses"] = [
            self.filter_definition["structure_thicknesses"][el]
            for el in self.layer_order
        ]
        temp_json["thickness_opt_allowed"] = [
            self.filter_definition["thickness_opt_allowed"][el]
            for el in self.layer_order
        ]
        temp_json["bounds"] = [
            self.filter_definition["bounds"][el] for el in self.layer_order
        ]
        temp_json["layer_switch_allowed"] = [
            self.filter_definition["layer_switch_allowed"][el]
            for el in self.layer_order
        ]
        temp_json["incoherent"] = [
            self.filter_definition["incoherent"][el] for el in self.layer_order
        ]

        return temp_json

    #####################################################
    ############# Filter Optimization Area ##############
    #####################################################

    def calculate_initial_merit(self):
        """
        Function to only calculate the current merit that is accessible from
        outside
        """
        my_filter, lib = self.create_filter_in_cpp()

        lib.initialise_optimization(
            my_filter,
            self.target_wavelength,
            int(np.size(self.target_wavelength)),
        )
        merit = lib.calculate_merit(
            my_filter,
            (ctypes.c_char_p * len(self.target_type.tolist()))(
                *[s.encode("utf-8") for s in self.target_type.tolist()]
            ),
            (ctypes.c_char_p * len(self.target_polarization.tolist()))(
                *[s.encode("utf-8") for s in self.target_polarization.tolist()]
            ),
            self.target_value,
            self.target_wavelength,
            self.target_polar_angle,
            self.target_azimuthal_angle,
            self.target_weights,
            (ctypes.c_char_p * len(self.target_condition.tolist()))(
                *[s.encode("utf-8") for s in self.target_condition.tolist()]
            ),
            self.target_tolerance,
            (ctypes.c_char_p * len(self.target_arithmetic.tolist()))(
                *[s.encode("utf-8") for s in self.target_arithmetic.tolist()]
            ),
            int(np.size(self.target_value)),
            # self.filter_definition["core_selection"],
        )

        return merit

    def perform_optimisation(
        self, opt_methods, redis_key=None, additional_parameters=[]
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
        opt_methods (str): The type of optimization to perform. This can
        be "dual_annealing", "differential_evolution", "basinhopping", "brute",
        "shgo", or "minimize".
        save_optimized_to_file (bool, optional): Whether to save the optimized
        values to a file. Defaults to True.
        stop_flag (function, optional): A function that returns True if the
        optimization should be stopped, False otherwise. Defaults to None.

        Returns:
        numpy array: The optimized features.
        """
        self.redis_key = redis_key

        if self.redis_key != None:
            # Setup Redis connection
            redis_url = os.environ.get(
                "REDIS_URL", "redis://localhost:6379"
            )  # Default to localhost if not set
            self.redis_conn = Redis.from_url(redis_url)

        log("loading filter stack...")

        my_filter, lib = self.create_filter_in_cpp()
        # self.job_id = job_id

        log("running optimisation...")

        iterator = 0

        for optimization_method in opt_methods:

            log("running " + optimization_method + "...")

            lib.initialise_optimization(
                my_filter,
                self.target_wavelength,
                int(np.size(self.target_wavelength)),
            )

            # Reset some variables
            self.stop_flag = False
            self.initial_merit = 0
            self.iteration_no = 0
            self.callback_call = 0
            self.optimization_method = optimization_method

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

            # Obtain the index of the layers that are allowed to switch
            # which are all where the layer_switch_allowed variable was set
            # to true and are coherent layers.
            switched_layers = np.where(
                np.logical_and(
                    self.filter_definition["layer_switch_allowed"],
                    ~np.array(self.filter_definition["incoherent"]),
                )
            )[0].tolist()

            if len(switched_layers) > 0:
                # If layer switching is allowed add the additional parameter that
                # decides the layer positions

                # Find the current position of these layers in the stack
                initial_positions = [
                    self.layer_order.tolist().index(element)
                    for element in switched_layers
                    if element in self.layer_order.tolist()
                ]

                x_initial = x_initial + initial_positions

                # Set the first layer positions to the initial positions
                self.previous_layer_positions = initial_positions

                for i in range(np.size(switched_layers)):
                    bounds.append(
                        (
                            -0.5,
                            np.size(self.initial_structure_thicknesses) + 0.5,
                        )
                    )

            ret = 0

            # Make sure the initial values are within the bounds otherwise the
            # results will be bad an hard to debug
            for x, bound in zip(x_initial, bounds):
                if not (bound[0] <= x <= bound[1]):
                    raise ValueError(
                        f"Initial value {x} is not within the bounds {bound}"
                    )

            partial_merit_function = partial(
                self.merit_function, my_filter=my_filter, lib=lib
            )

            # With scipy we cannot do integer optimization
            if optimization_method == "dual_annealing":
                ret = dual_annealing(
                    partial_merit_function,
                    bounds=bounds,
                    callback=self.scipy_callback,
                    x0=x_initial,
                    maxiter=10000,
                    minimizer_kwargs={"callback": self.scipy_callback},
                )
            elif optimization_method == "differential_evolution":
                ret = differential_evolution(
                    partial_merit_function,
                    bounds=bounds,
                    x0=x_initial,
                    maxiter=100000,
                    callback=self.scipy_callback,
                )
            elif optimization_method == "basinhopping":
                # This algorithm does
                # 1. a random perturbation of the features
                # 2. then a scipy.minimize
                # 3. accepts of rejects the new optimum value
                ret = basinhopping(
                    partial_merit_function,
                    x0=x_initial,
                    callback=self.scipy_callback,
                    minimizer_kwargs={
                        "method": "Nelder-Mead",
                        "bounds": bounds,
                        "callback": self.scipy_callback,
                    },
                )

            elif optimization_method == "brute":
                # Brute force optimization: only sensible for a low number of
                # features (e.g., 3). Depending on Ns, the number of points to
                # try is established  (Ns^(#features) = number of iterations)
                ret = brute(
                    partial_merit_function,
                    ranges=bounds,
                    Ns=2,
                    # maxiter = 50000,
                )
            elif optimization_method == "shgo":
                # Doesn't really start
                ret = shgo(
                    partial_merit_function,
                    bounds=bounds,
                    # maxiter = 50000,
                )
            elif optimization_method == "LM":
                ret = gradient_descent(
                    partial_merit_function,
                    x0=x_initial,
                    bounds=bounds,
                    callback=self.scipy_callback,
                )
            elif optimization_method == "TNC":
                # Truncated newton method
                ret = minimize(
                    partial_merit_function,
                    x0=x_initial,
                    bounds=bounds,
                    method="TNC",  # 41570 after 10000 iterations
                    jac=lambda x: gradient(self.merit_function, x),
                    callback=self.scipy_callback,
                )
            elif optimization_method == "Nelder-Mead":
                # Nelder-Mead method
                ret = minimize(
                    partial_merit_function,
                    x0=x_initial,
                    bounds=bounds,
                    method="Nelder-Mead",
                    callback=self.scipy_callback,
                )
            elif optimization_method == "particle swarm":
                # Particle swarm optimization
                ret = particle_swarm(
                    partial_merit_function,
                    bounds=bounds,
                    n_particles=additional_parameters[iterator]["particles"],
                    c1=additional_parameters[iterator]["c1"],
                    c2=additional_parameters[iterator]["c2"],
                    w=additional_parameters[iterator]["w"],
                    n_iter=additional_parameters[iterator]["iterations"],
                    callback=self.scipy_callback,
                )
            elif optimization_method == "grid search particle swarm hyperparameter":
                # Particle swarm optimization
                grid_search(
                    partial_merit_function,
                    bounds=bounds,
                    n_particles=additional_parameters[iterator]["particles"],
                    c1=additional_parameters[iterator]["c1"],
                    c2=additional_parameters[iterator]["c2"],
                    w=additional_parameters[iterator]["w"],
                    n_iter=additional_parameters[iterator]["iterations"],
                    callback=self.scipy_callback,
                )
            else:
                raise ValueError(
                    f"Optimization method {optimization_method} not implemented"
                )

            """
            # Wrap the arguments to pass them to the Process
            def worker(
                optimization_method, partial_merit_function, x_initial, bounds, callback
            ):
                optimize_in_process(
                    optimization_method,
                    partial_merit_function,
                    x_initial,
                    bounds,
                    callback,
                )

            # Create and start a separate process for the optimization task
            p = Process(
                target=worker,
                args=(
                    optimization_method,
                    partial_merit_function,
                    x_initial,
                    bounds,
                    self.scipy_callback,
                ),
            )
            p.start()
            p.join()  # Wait for the optimization to complete if necessary
            """

            # The result from the optimization is not necessarily the global
            # result (as e.g. gradients are also calculated that could have a
            # lower merit)
            ret.x = self.optimum_x
            ret.fun_best = self.optimum_merit

            thicknesses, self.layer_order = (
                self.extract_thickness_and_position_from_features(self.optimum_x)
            )

            # Change the filter stack to the optimized values (in C++ and in python)
            self.filter_definition["structure_thicknesses"] = thicknesses
            lib.change_material_thickness(
                my_filter, thicknesses, int(np.size(thicknesses))
            )
            self.filter_definition["structure_materials"] = [
                self.initial_structure_materials[el] for el in self.layer_order
            ]
            lib.change_material_order(
                my_filter, self.layer_order, int(np.size(self.layer_order))
            )

            lib.get_material_order(my_filter)
            lib.get_thicknesses(my_filter)

            # The stop flag is also be used as an "optimization done" indicator
            self.stop_flag = True

            self.log_func("Optimization time: " + str(time.time() - start_time) + "s")
            self.log_func(
                "Optimized thicknesses: "
                + str([thicknesses[el] for el in self.layer_order])
            )
            self.log_func(
                "Optimized layer order: "
                + str(self.filter_definition["structure_materials"])
            )
            self.log_func("Optimized merit value: " + str(self.optimum_merit))
            self.log_func("Number of function evaluations: " + str(ret.nfev))
            iterator += 1

        """
        if key != None:
            # Close the Redis connection explicitly
            self.redis_conn.close()
        """

        return ret.x

    def merit_function(self, features, my_filter, lib):
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

        # Extract thickness and layer order from features
        thicknesses, self.layer_order = (
            self.extract_thickness_and_position_from_features(features)
        )

        # Change material thickness
        if np.any(self.filter_definition["thickness_opt_allowed"]):
            self.filter_definition["structure_thicknesses"] = thicknesses
            lib.change_material_thickness(
                my_filter, thicknesses, int(np.size(thicknesses))
            )

        # Change layer switch allowed

        if np.any(
            np.logical_and(
                self.filter_definition["layer_switch_allowed"],
                ~np.array(self.filter_definition["incoherent"]),
            )
        ):
            self.filter_definition["structure_materials"] = [
                self.initial_structure_materials[el] for el in self.layer_order
            ]
            lib.change_material_order(
                my_filter, self.layer_order, int(np.size(self.layer_order))
            )

        # Calculate merit
        merit = lib.calculate_merit(
            my_filter,
            (ctypes.c_char_p * len(self.target_type.tolist()))(
                *[s.encode("utf-8") for s in self.target_type.tolist()]
            ),
            (ctypes.c_char_p * len(self.target_polarization.tolist()))(
                *[s.encode("utf-8") for s in self.target_polarization.tolist()]
            ),
            self.target_value,
            self.target_wavelength,
            self.target_polar_angle,
            self.target_azimuthal_angle,
            self.target_weights,
            (ctypes.c_char_p * len(self.target_condition.tolist()))(
                *[s.encode("utf-8") for s in self.target_condition.tolist()]
            ),
            self.target_tolerance,
            (ctypes.c_char_p * len(self.target_arithmetic.tolist()))(
                *[s.encode("utf-8") for s in self.target_arithmetic.tolist()]
            ),
            int(np.size(self.target_value)),
            # self.filter_definition["core_selection"],
        )

        # If merit is zero, the optimization has reached the target
        if merit == 0:

            # Callback once to save the obtained features
            if self.iteration_no == 0:
                raise Exception("Initial features are optimal.")

            if self.first_zero:
                self.first_zero = False
                self.log_func("Optimization has reached merit 0")
                self.callback(features, merit)
                self.stop_flag = True

            return 0

        else:

            # Set initial merit and normalize to it
            if self.initial_merit == 0:
                self.initial_merit = merit
                self.optimum_merit = merit

            ## callback function
            self.callback(features, merit)

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
        intermediate_result = {
            "step": self.iteration_no,
            "merit": f,
        }
        if self.redis_key != None:
            self.redis_conn.set(
                f"current:{self.redis_key}",
                json.dumps(intermediate_result),
            )

        # Save the current best optimisation values to file
        if f < self.optimum_merit or f == 0:

            # in webapp, current_structure will be linked to the session_id
            if self.redis_key != None:
                intermediate_result = {
                    "step": self.iteration_no,
                    "merit": f,
                    "current_structure": self.return_current_design_as_json(),
                }
                self.redis_conn.set(
                    f"current_best:{self.redis_key}",
                    json.dumps(intermediate_result),
                )
            else:
                file_name = self.current_structure.split("/")[-1]
                temp_json = self.return_current_design_as_json()

                temp_dir = os.path.join(
                    os.path.dirname(os.path.realpath(__file__)), "temp"
                )

                # Create the directory if it does not exist
                os.makedirs(temp_dir, exist_ok=True)

                temp_path = os.path.join(temp_dir, f"{file_name}.json")

                with open(temp_path, "w") as file:
                    json.dump(temp_json, file)
            """
            if database_entry is not None:
                # Assuming database_entry is a single SQLAlchemy model instance
                database_entry.current_structure = self.current_structure
                database_entry.current_merit = f
                database_entry.current_iteration = self.iteration_no

                db.session.commit()  # Commit the changes to the database
            """
            """
            if self.current_structure != "current_structure":
                temp_path = os.path.join(
                    os.path.dirname(os.path.realpath(__file__)),
                    "temp",
                    f"{self.current_structure.split('/')[-1]}.json",
                )
                with open(f"{self.current_structure}.pkl", "wb") as file_pickled:
                    data_to_pickle = (
                        self,
                        self.job_id,
                    )
                    pickle.dump(data_to_pickle, file_pickled)
            """

            self.optimum_merit = f
            self.optimum_iteration = self.iteration_no
            self.optimum_number += 1
            self.optimum_x = x

        self.callback_call += 1

        current_time = time.time()

        # Number of seconds to wait for next logging is defined here at 2 seconds
        if current_time - self.last_log_time >= 2:
            self.log_func(
                f"merit | call #{str(self.iteration_no)} : {round(f / self.initial_merit, 9)} {round(f, 2)}"
            )
            self.last_log_time = current_time
            self.last_merit = f
            if self.last_optimum_number < self.optimum_number:
                self.log_func(
                    f"New optimum on call #{self.optimum_iteration} : {round(self.optimum_merit / self.initial_merit, 9)} {round(self.optimum_merit, 2)}"
                )
                self.log_design_func()
                self.last_optimum_number = self.optimum_number

        # Implement cancelation depending on the job_cancel_flag in the redis database
        if self.redis_key != None:
            cancel_flag = self.redis_conn.get(f"job_cancel_flag:{self.redis_key}")

            if (
                cancel_flag is not None
                and cancel_flag.decode("utf-8").lower() == "true"
            ):
                # Perform any necessary cleanup here
                self.stop_flag = True
                print("CANCELED!")

        # If the merit function is close to zero (within the tolerance) or the
        # stop flag is true, stop the optimization.
        # if self.stop_flag or math.isclose(f, 0, abs_tol=1e-6):
        # return True
        # else:
        # return False

    def scipy_callback(self, dummy1, dummy2=None, dummy3=None):
        """
        This is used to stop the scipy optimize function
        """
        if self.stop_flag:
            if self.optimization_method == "TNC":
                raise StopIteration
            # elif self.optimization_method == "NCG":
            # raise StopIteration
            elif self.optimization_method == "Nelder-Mead":
                raise StopIteration
            elif self.optimization_method == "LM":
                return True
            elif self.optimization_method == "dual_annealing":
                return True
            elif self.optimization_method == "basinhopping":
                return True
            elif self.optimization_method == "differential_evolution":
                return True
            elif self.optimization_method == "particle swarm":
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
        layer_switch_allowed_without_incoherent = np.logical_and(
            self.filter_definition["layer_switch_allowed"],
            ~np.array(self.filter_definition["incoherent"]),
        )

        if np.any(self.filter_definition["thickness_opt_allowed"]) and not np.any(
            layer_switch_allowed_without_incoherent
        ):
            # All features are thicknesses
            thicknesses = np.copy(self.initial_structure_thicknesses)
            thicknesses[np.where(self.filter_definition["thickness_opt_allowed"])] = (
                features
            )
            thicknesses = thicknesses.astype(np.float64)
            layer_order = np.arange(0, len(self.initial_structure_thicknesses))
        elif not np.any(self.filter_definition["thickness_opt_allowed"]) and np.any(
            layer_switch_allowed_without_incoherent
        ):
            # All features are layer positions
            thicknesses = np.copy(self.initial_structure_thicknesses)

            # Transform the layer positions to a full stack order
            layer_order = self.convert_layer_positions_to_stack_order(features)

            # layer_order = self.allowed_permutations[
            # self.clamp(features[-1], 0, len(self.allowed_permutations) - 1)
            # ].astype(np.int32)
        elif np.any(self.filter_definition["thickness_opt_allowed"]) and np.any(
            layer_switch_allowed_without_incoherent
        ):
            # Some features are thicknesses and some are layer positions
            thicknesses = np.copy(self.initial_structure_thicknesses)
            thicknesses[np.where(self.filter_definition["thickness_opt_allowed"])] = (
                features[: -1 * np.sum(layer_switch_allowed_without_incoherent)]
            )
            thicknesses = thicknesses.astype(np.float64)

            # The hard part is to find the correct way of generating a sensible
            # structure from the feature numbers
            temp_layer_positions = features[
                -1 * np.sum(layer_switch_allowed_without_incoherent) :
            ]

            # Transform the layer positions to a full stack order
            layer_order = self.convert_layer_positions_to_stack_order(
                temp_layer_positions
            )

        else:
            raise ValueError

        # Clip the thicknesses to the bounds (but only if an optimization was done in the first place)
        return np.array(
            [
                (
                    np.clip(np.nan_to_num(np.round(t, 1), nan=0.0), b[0], b[1])
                    if self.filter_definition["thickness_opt_allowed"][i]
                    else t
                )
                for i, (t, b) in enumerate(zip(thicknesses, self.bounds))
            ],
            dtype=np.float64,
        ), layer_order.astype(np.int32)

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
        temp_layer_positions = np.round(temp_layer_positions, 0)

        # Now find out if some incoherent layers were set to contain
        containment_layers = np.where(
            np.logical_and(
                self.filter_definition["layer_switch_allowed"],
                self.filter_definition["incoherent"],
            )
        )[0].tolist()

        # Now clamp the integers to either their containment or to the available
        # number of layers
        if len(containment_layers) == 0:
            temp_layer_positions = np.clip(
                temp_layer_positions,
                0,
                len(self.initial_structure_thicknesses),
            ).astype(np.int32)
        else:
            # Find out which of the layers belong to which containment zone
            # bracketed by the containment layer and the indexes zero and the
            # maximum index
            initial_switchable_layer_positions = np.where(
                np.logical_and(
                    self.filter_definition["layer_switch_allowed"],
                    ~np.array(self.filter_definition["incoherent"]),
                )
            )[0]

            for i in range(len(containment_layers)):
                # If above the containment layer, clamp to higher than its position
                temp_layer_positions[
                    np.where(initial_switchable_layer_positions > containment_layers[i])
                ] = np.clip(
                    temp_layer_positions[
                        np.where(
                            initial_switchable_layer_positions > containment_layers[i]
                        )
                    ],
                    containment_layers[i] + 1,
                    len(self.initial_structure_thicknesses),
                ).astype(
                    np.int32
                )

                # If below containment layer, clamp to lower than its position
                temp_layer_positions[
                    np.where(initial_switchable_layer_positions < containment_layers[i])
                ] = np.clip(
                    temp_layer_positions[
                        np.where(
                            initial_switchable_layer_positions < containment_layers[i]
                        )
                    ],
                    0,
                    containment_layers[i] - 1,
                ).astype(
                    np.int32
                )

        # If feature 1 has a certain number, the second layer cannot have
        # the same (one position can only be occupied by one layer)
        unique_layer_positions = self.ensure_unique(temp_layer_positions)

        # Set previous layer positions to be the current ones
        self.previous_layer_positions = unique_layer_positions

        # Start with the initial order of layers
        layer_order = list(range(len(self.initial_structure_thicknesses)))

        # Iterate over each layer
        j = 0
        for i in range(len(self.filter_definition["layer_switch_allowed"])):
            if (
                self.filter_definition["layer_switch_allowed"][i]
                and not self.filter_definition["incoherent"][i]
            ):
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
        return arr.astype(np.int32)

    #####################################################
    ########## Plotting and calculation area ############
    #####################################################

    def read_results_from_cpp(self, result_string, wavelengths, polar_angles):
        """
        Reads the results from the C++ calculation and returns them as a
        pandas DataFrame.
        """
        # List of pandas dataframes with one data frame per azimuthal angle
        stored_values = []

        # Split the string by "=" to obtain data by azimuthal angle
        azim_angle_split_arrays = result_string.decode("utf-8").split("=")

        # Iterate over azimuthal angles
        for azim_angle_split_array in azim_angle_split_arrays:
            stored_value = pd.DataFrame(
                columns=np.ascontiguousarray(polar_angles).astype(np.float64),
                index=wavelengths,
            )

            # The polar angle columns are separated by "+", so split for them
            theta_angle_split_arrays = azim_angle_split_array.split("+")

            # Iterate over all polar angles and split again (by "--" this time
            # to obtain data for each wavelength
            for n, theta_angle_split_array in enumerate(theta_angle_split_arrays):
                stored_value.loc[:, polar_angles[n]] = np.array(
                    theta_angle_split_array.split("--"),
                    dtype=np.float64,
                )
            stored_values.append(stored_value)

        return stored_values

    """
    # Not used at the moment!
    def calculate_one_angle(
        self,
        minimum_wavelength,
        maximum_wavelength,
        wavelength_step,
        target_type,
        polarization,
        polar_angle,
        azimuthal_angle,
        # is_general_core,
    ):
        # Function to only get the results for one particular angle given a wavelength range

        my_filter, lib = self.create_filter_in_cpp()

        wavelengths = np.arange(
            minimum_wavelength, maximum_wavelength + 1, wavelength_step
        )

        if polarization == "s":
            polarization = 1.0
        elif polarization == "p":
            polarization = 0.0
        elif polarization == "":
            polarization = 0.5
        elif isinstance(polarization, (int, float)) and 0 <= polarization <= 1:
            # Lastly check if polarization is numeric and between [0, 1]
            pass
        else:
            raise ValueError("Polarization must be either '', 's', 'p', or given as s-polarization as a number between 0 and 1")

        result_string = lib.calculate_reflection_transmission_absorption_para(
            my_filter,
            target_type.encode("utf-8"),
            polarization.astype(float),
            np.ascontiguousarray(wavelengths).astype(np.float64),
            int(np.size(wavelengths)),
            np.ascontiguousarray(polar_angle).astype(np.float64),
            int(np.size(polar_angle)),
            np.ascontiguousarray(azimuthal_angle).astype(np.float64),
            int(np.size(azimuthal_angle)),
            # is_general_core,
        )

        return pd.Series(
            np.array(result_string.decode("utf-8").split("--"), dtype=np.float64),
            index=np.ascontiguousarray(wavelengths).astype(np.float64),
        )
    """

    def calculate_ar_data(
        self,
        wavelength=None,
        polar_angles=None,
        azimuthal_angles=None,
        target_type=None,
        polarization=None,
        save_figure=False,
        save_data=False,
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

        my_filter, lib = self.create_filter_in_cpp()

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

        if polarization == "s":
            polarization = 1.0
        elif polarization == "p":
            polarization = 0.0
        elif polarization == "":
            polarization = 0.5
        elif isinstance(polarization, (int, float)) and 0 <= polarization <= 1:
            # Lastly check if polarization is numeric and between [0, 1]
            pass
        else:
            raise ValueError(
                "Polarization must be either '', 's', 'p', or given as s-polarization as a number between 0 and 1"
            )

        initial_time = time.time()

        result_string = lib.calculate_reflection_transmission_absorption_para(
            my_filter,
            target_type.encode("utf-8"),
            float(polarization),
            np.ascontiguousarray(wavelength).astype(np.float64),
            int(np.size(wavelength)),
            np.ascontiguousarray(polar_angles).astype(np.float64),
            int(np.size(polar_angles)),
            np.ascontiguousarray(azimuthal_angles).astype(np.float64),
            int(np.size(azimuthal_angles)),
        )

        self.stored_data = self.read_results_from_cpp(
            result_string, wavelength, polar_angles
        )

        if save_data:
            header_lines = []

            # Save header lines indicating what the simulation represents
            temp_path = os.path.join(
                os.path.dirname(os.path.realpath(__file__)), "./temp/value.csv"
            )
            # temp_path = os.path.join(os.getcwd(), "temp", "value.csv")
            with open(temp_path, "w") as the_file:
                the_file.write("\n".join(header_lines))

            # Save actual data by appending
            self.stored_data[0].to_csv(temp_path, sep=",", header=True, mode="a")
        if save_figure:
            # Print time elapsed for the generation of the reflectivity matrix
            log(time.time() - initial_time)
            # Plotting - right now, on the first azimuthal angle
            plt.close()
            X, Y = np.meshgrid(
                self.stored_data[0].columns.astype(float),
                self.stored_data[0].index.astype(float),
            )
            # Prepare data for 3D plot where each column contains the same data for
            # the different angles
            Z = self.stored_data[0].to_numpy(float)
            plt.pcolormesh(X, Y, Z, shading="auto")
            # Add a colorbar
            plt.colorbar(label="Intensity (a.u.)")
            # Visuals
            plt.xlabel("Polar Angle (°)")
            plt.ylabel("Wavelength (nm)")
            # Save the figure before showing it
            # plt.savefig(f"{phi}-plot.png", format="png", dpi=300)
            temp_path = os.path.join(
                os.path.dirname(os.path.realpath(__file__)), "./temp/plot.png"
            )
            # temp_path = os.path.join(os.path.dirname(os.getcwd()), "temp", "plot.png")
            plt.savefig(temp_path, format="png", dpi=300)
            plt.show()
            # Save X, Y, Z to csv files
