from FilterStack import FilterStack
import json

## typical call: python main.py dual_annealing ../examples/AlPcCl_bandpass_840nm.json
## default values are current_structure.json from last optimisation and dual_annealing

"""
# Create the parser
parser = argparse.ArgumentParser()

# Add the arguments
parser.add_argument(
    "optimisation_method",
    type=str,
    help="The optimisation method",
    default=["LM", "dual_annealing"],
    nargs="?",
)
parser.add_argument(
    "file",
    type=str,
    help="The input file path",
    default="../examples/demo_test.json",
    nargs="?",
)

# Parse the arguments
args = parser.parse_args()
"""
json_path = "../examples/demo_test.json"
optimisation_method = ["particle swarm"]
additional_parameters = [
    {"particles": 30, "c1": 1.5, "c2": 1.5, "w": 2.5, "iterations": 10000}
]

# Read the JSON file
with open(json_path, "r") as file:
    filter_json = json.load(file)

# Create filter stack
filter_stack = FilterStack(my_filter_dict=filter_json)

features = filter_stack.perform_optimisation(
    optimisation_method, additional_parameters=additional_parameters
)
