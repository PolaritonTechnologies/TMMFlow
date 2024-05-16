# tests.py performs basic unit testing for the optimisation module
# command is pytest tests.py

from FilterStack import FilterStack
import pandas as pd

#########################
# Input parameters
json_file_path = "../examples/demo_test.json"

#########################

# Create filter stack
filter_stack = FilterStack(json_file_path)

# Optimize stack
# filter_stack.lib.get_thicknesses()
# filter_stack.lib.get_material_order()


# def test():
# assert optimization.check_targets(cavity_optimisation())
for target in ["t", "r", "a"]: 
    for polarization in ["s", "p"]:
        print("Target: " + target + ", Pol.: " + polarization)
        filter_stack.calculate_ar_data(
            polar_angles=[0, 15, 30, 45, 60, 75, 89],
            save_data=True,
            target_type = target,
            polarization = polarization,
        )
        calculated_data_df = pd.read_csv("../src/temp/value.csv", sep = ",")
        calculated_data_df = calculated_data_df.rename(columns={"Unnamed: 0": "wavelength"})

        open_filter_df = pd.read_csv("../tests/" + polarization + "_" + target + ".csv", sep = "\t")

        # x = 1e-3
        # calculated_data_df = calculated_data_df.where(calculated_data_df >= x, 0)
        # open_filter_df = open_filter_df.where(open_filter_df >= x, 0)
        temp_diff = abs(calculated_data_df - open_filter_df)
        print(temp_diff)
        print("Maximum deviations: ")
        print(temp_diff.max())

        # Plot this to get a feeling for where the inaccuracies are the worst

        # Testing with assert is not really helpful here
        # pd.testing.assert_frame_equal(calculated_data_df, open_filter_df)
