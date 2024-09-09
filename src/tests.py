# tests.py performs basic unit testing for the optimisation module
# command is pytest tests.py

from FilterStack import FilterStack
import pandas as pd
import numpy as np
import json
import os

#########################
# Convert data exported from openFilters for comparison to our standard
# dataframe format


def convert_open_filter_datastructure(file_path, columns):
    """
    Convert the weird open filter file format to something comparable to our
    pandas dataframe
    """
    # Read the file line by line, skipping lines that start with whitespace
    # followed by 'wavelength', 'Reflection', 'Transmission' or 'Absorption'
    with open(file_path, "r") as f:
        lines = [
            line
            for line in f
            if not line.lstrip().startswith(
                ("wavelength", "Reflection", "Transmission", "Absorption")
            )
        ]

    # Convert the list of lines to a single string
    data = "\n".join(lines)

    # Use StringIO to read the string as a file
    from io import StringIO

    data_io = StringIO(data)

    # Read the CSV data
    df = pd.read_csv(data_io, sep="\s+", names=["wavelength", "R"])

    # Reorder the dataframe so that it contains
    number_of_angles = int(np.size(columns) - 1)
    unique_wavelength = np.unique(df.to_numpy().T[0])
    intensity_data = (
        df.to_numpy()
        .T[1]
        .reshape(int(number_of_angles), int(df.to_numpy().T[1].size / number_of_angles))
    )
    full_frame_data = np.concatenate(
        [unique_wavelength.reshape(1, np.size(unique_wavelength)), intensity_data],
        axis=0,
    )
    df_to_write = pd.DataFrame(full_frame_data.T, columns=columns)

    # Write dataframe to file
    df_to_write.to_csv(
        ".".join(np.append(file_path.split(".")[:-1], "_conv.csv")),
        sep="\t",
        index=False,
    )


# convert_open_filter_datastructure("../tests/a_p_flat_top.csv", ['wavelength', 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0])
# convert_open_filter_datastructure("../tests/a_s_flat_top.csv", ['wavelength', 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0])
# convert_open_filter_datastructure("../tests/t_p_flat_top.csv", ['wavelength', 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0])
# convert_open_filter_datastructure("../tests/t_s_flat_top.csv", ['wavelength', 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0])
# convert_open_filter_datastructure("../tests/r_p_flat_top.csv", ['wavelength', 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0])
# convert_open_filter_datastructure("../tests/r_s_flat_top.csv", ['wavelength', 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0])

#########################
# Input parameters
# json_file_path = os.path.join(
# os.path.dirname(os.getcwd()), "examples", "demo_test_with_backside.json"
# )


def execute_tests(json_path):

    json_file_path = os.path.join(os.path.dirname(os.getcwd()), "examples", json_path)

    #########################

    with open(json_file_path, "r") as file:
        filter_json = json.load(file)

    # Create filter stack
    filter_stack = FilterStack(my_filter_dict=filter_json)

    # def test():
    # assert optimization.check_targets(cavity_optimisation())
    for target in ["t", "r"]:
        for polarization in ["s", "p"]:
            print("Target: " + target + ", Pol.: " + polarization)
            filter_stack.calculate_ar_data(
                polar_angles=[0, 15, 30, 45, 60, 75, 89],
                save_data=True,
                target_type=target,
                polarization=polarization,
            )
            calculated_data_df = pd.read_csv(
                os.path.join(os.getcwd(), "temp", "value.csv"), sep=","
            )
            calculated_data_df = calculated_data_df.rename(
                columns={"Unnamed: 0": "wavelength"}
            )

            # To find the corresponding test file, one needs to remove json from the filename
            cleaned_file_name = json_path.replace(".json", "")

            open_filter_df = pd.read_csv(
                os.path.join(os.path.dirname(os.getcwd()), "src/tests")
                + "/"
                + cleaned_file_name
                + "_"
                + polarization
                + "_"
                + target
                + ".csv",
                sep="\t",
            )

            # complete_ease_df = pd.read_csv("../tests/" + polarization + "_" + target + "_CE.csv", sep = "\t")

            # x = 1e-3
            # calculated_data_df = calculated_data_df.where(calculated_data_df >= x, 0)
            # open_filter_df = open_filter_df.where(open_filter_df >= x, 0)
            temp_diff = abs(calculated_data_df - open_filter_df)
            # temp_diff_ce_of = abs(complete_ease_df - open_filter_df)
            # print(temp_diff)
            # print(temp_diff_ce_of)
            print("Maximum deviations: ")
            print(temp_diff.max())
            # print(temp_diff_ce_of.max())

            # Plot this to get a feeling for where the inaccuracies are the worst

            # Testing with assert is not really helpful here
            # pd.testing.assert_frame_equal(calculated_data_df, open_filter_df)


tests = ["demo_test.json", "demo_test_with_backside.json"]

for test in tests:
    execute_tests(test)
