import argparse
from FilterStack import FilterStack
import numpy as np

## typical call: python main.py dual_annealing ../examples/AlPcCl_bandpass_840nm.json
## default values are current_structure.json from last optimisation and dual_annealing

# Create the parser
parser = argparse.ArgumentParser()

# Add the arguments
parser.add_argument(
    "optimisation_method",
    type=str,
    help="The optimisation method",
    default="minimize",
    nargs="?",
)
parser.add_argument(
    "file",
    type=str,
    help="The input file path",
    default="../src/temp/current_structure.json",
    nargs="?",
)

# Parse the arguments
args = parser.parse_args()

# Create filter stack
filter_stack = FilterStack(args.file)

# Optimize stack
# filter_stack.lib.get_thicknesses()
# filter_stack.lib.get_material_order()

features = filter_stack.perform_optimisation(
    args.optimisation_method, save_optimized_to_file=True
)
# filter_stack.merit_function(np.append(filter_stack.filter_definition["structure_thicknesses"],[3,19,21]))
# print(filter_stack.merit_function(features) *  filter_stack.initial_merit)

# Calculate the data for
filter_stack.calculate_ar_data(
    save_figure=True,
    save_data=True,
)
