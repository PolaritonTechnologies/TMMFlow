import argparse
from FilterStack import FilterStack

## typical call: python main.py ../examples/AlPcCl_bandpass_840nm.json dual_annealing
## default values are current_structure.json from last optimisation and dual_annealing

# Create the parser
parser = argparse.ArgumentParser()

# Add the arguments
parser.add_argument(
    "--file",
    type=str,
    help="The input file path",
    default="./temp/current_structure.json",
)
parser.add_argument(
    "--optimisation_method",
    type=str,
    help="The optimisation method",
    default="dual_annealing",
)

# Parse the arguments
args = parser.parse_args()

# Create filter stack
filter_stack = FilterStack(args.file)

# Optimize stack
features = filter_stack.perform_optimisation(
    args.optimisation_method, save_optimized_to_file=True
)

# Calculate the data for
filter_stack.calculate_ar_data(
    save_figure=True,
    save_data=True,
)
