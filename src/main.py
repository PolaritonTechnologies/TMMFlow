from FilterStack import FilterStack

# Input file
filter_definition_json_file_path = "../examples/AlPcCl_bandpass_840nm.json"

# Create filter stack
filter_stack = FilterStack(filter_definition_json_file_path)

# Optimize stack
features = filter_stack.perform_optimisation(
    "dual_annealing", save_optimized_to_file=True
)

# Calculate the data for 
filter_stack.calculate_ar_data(
    save_figure=True,
    save_data=True,
)
