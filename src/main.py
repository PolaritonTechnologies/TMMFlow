import ctypes
import json
import numpy as np

from optim_module import OptimModule
from calculation_module import CalculationModule

from utility import translate_order_for_cpp, create_filter

#########################
# Input parameters
optimisation_order_file_python = "demo_test.json"
optimisation_order_file = translate_order_for_cpp(optimisation_order_file_python)
my_filter, lib = create_filter(optimisation_order_file)

#########################
# Optimization
print("running optimisation...")
optimization = OptimModule(optimisation_order_file, my_filter, lib)
features = optimization.perform_optimisation(
    "minimize", save_optimized_to_file=True, stop_flag=lambda: False
)
#########################

#########################

# calculation_order_file_python = "demo_test.json"
# calculation_order_file = translate_order_for_cpp(calculation_order_file_python)
# my_filter, lib = create_filter(calculation_order_file)

# Plotting - Can also implement core override here if needed
plotting = CalculationModule(my_filter, lib, core_selection="general", web=False)

with open(optimisation_order_file) as f:
    plot_order = json.load(f)

wavelength = np.arange(
    plot_order["wavelengthMin"],
    plot_order["wavelengthMax"] + 1,
    plot_order["wavelengthStep"],
)
polar_angles = np.arange(
    plot_order["polarAngleMin"],
    plot_order["polarAngleMax"] + 1,
    plot_order["polarAngleStep"],
)
azimuthal_angles = np.arange(
    plot_order["azimAngleMin"],
    plot_order["azimAngleMax"] + 1,
    plot_order["azimAngleStep"],
)

print("plotting results...")
plotting.calculate_ar_data(
    wavelength,
    polar_angles,
    azimuthal_angles,
    plot_order["calculation_type"],
    plot_order["polarization"],
    save_figure=True,
    save_data=True,
)
#########################
