# tests.py performs basic unit testing for the optimisation module
# command is pytest tests.py

from FilterStack import FilterStack
import pandas as pd
import numpy as np
import json
import os
import time


#########################
# Convert data exported from openFilters for comparison to our standard
# dataframe format
#########################
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


def execute_tests(test_dic, verbosity="loud"):

    json_file_path = os.path.join(
        os.path.dirname(os.getcwd()), "examples", test_dic["json_path"]
    )

    polarization_to_float = {"s": 1, "p": 0, "": 0.5}

    #########################

    if verbosity == "quiet":
        # redirect print to nothing
        import sys

        sys.stdout = open(os.devnull, "w")

    with open(json_file_path, "r") as file:
        filter_json = json.load(file)

    # Create filter stack
    filter_stack = FilterStack(my_filter_dict=filter_json)

    # Store the maximum deviations to return them later
    max_deviations_test = []

    for azimuthal_angle in test_dic["azimuthal_angles"]:

        for target in test_dic["calculation_types"]:

            for polarization in test_dic["polarizations"]:

                print(
                    "Target: "
                    + target
                    + ", Pol.: "
                    + polarization
                    + ", Az.: "
                    + str(azimuthal_angle)
                )

                filter_stack.calculate_ar_data(
                    polar_angles=test_dic["polar_angles"],
                    save_data=True,
                    target_type=target,
                    polarization=polarization_to_float[polarization],
                )

                calculated_data_df = pd.read_csv(
                    os.path.join(os.getcwd(), "temp", "value.csv"), sep=","
                )

                calculated_data_df = calculated_data_df.rename(
                    columns={"Unnamed: 0": "wavelength"}
                )

                # To find the corresponding test file, one needs to remove json from the filename
                cleaned_file_name = test_dic["json_path"].replace(".json", "")

                if polarization != "":

                    if test_dic["azimuthal_angles"] == [0]:

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

                    else:

                        open_filter_df = pd.read_csv(
                            os.path.join(os.path.dirname(os.getcwd()), "src/tests")
                            + "/"
                            + cleaned_file_name
                            + "_"
                            + str(azimuthal_angle)
                            + "deg"
                            + "_"
                            + polarization
                            + "_"
                            + target
                            + ".csv",
                            sep="\t",
                        )

                else:

                    open_filter_df_s = pd.read_csv(
                        os.path.join(os.path.dirname(os.getcwd()), "src/tests")
                        + "/"
                        + cleaned_file_name
                        + "_"
                        + "s"
                        + "_"
                        + target
                        + ".csv",
                        sep="\t",
                    )

                    open_filter_df_p = pd.read_csv(
                        os.path.join(os.path.dirname(os.getcwd()), "src/tests")
                        + "/"
                        + cleaned_file_name
                        + "_"
                        + "p"
                        + "_"
                        + target
                        + ".csv",
                        sep="\t",
                    )

                    # Combine the two dataframes with doing a sum of the columns with a factor 0.5 for each term
                    open_filter_df = open_filter_df_s.add(
                        open_filter_df_p, fill_value=0
                    )
                    open_filter_df = open_filter_df * 0.5

                temp_diff = abs(calculated_data_df - open_filter_df)
                print("Maximum deviations: ")
                print(temp_diff.max())
                # among the values of temp_diff.max() we need to find the maximum value
                max_deviations_test.append(temp_diff.max().max())

        if verbosity == "quiet":
            # redirect print back to stdout
            sys.stdout = sys.__stdout__

        return max_deviations_test


test_dic_demo_test = {
    "json_path": "demo_test.json",
    "polar_angles": [0, 15, 30, 45, 60, 75, 89],
    "azimuthal_angles": [0],
    "calculation_types": ["t", "r"],
    "polarizations": ["s", "p", ""],
}

test_dic_demo_test_with_backside = {
    "json_path": "demo_test_with_backside.json",
    "polar_angles": [0, 15, 30, 45, 60, 75, 89],
    "azimuthal_angles": [0],
    "calculation_types": ["t", "r"],
    "polarizations": ["s", "p", ""],
}

test_dic_pfo_aligned = {
    "json_path": "demo_test_pfo_aligned.json",
    "polar_angles": [0, 10, 20, 30, 40, 50, 60, 70],
    "azimuthal_angles": [0, 15, 30, 45, 60, 75, 90],
    "calculation_types": ["r"],
    "polarizations": ["s", "p"],
}

# test_dics = [test_dic_demo_test, test_dic_demo_test_with_backside, test_dic_pfo_aligned]
test_dics = [test_dic_demo_test, test_dic_demo_test_with_backside]
# test_dics = [test_dic_pfo_aligned]

max_deviations = []
verbosity = "loud"

if verbosity == "loud":
    print("Running tests...")
else:
    print("Running tests quietly...")

# start the timer
start = time.time()

for test in test_dics:
    max_deviations.append(execute_tests(test, verbosity=verbosity))

# find the maximum among the maximum deviations
max_dev_it = max(max_deviations)
# find the max inside the max dev list
max_dev = max(max_dev_it)

max_acceptable = 0.003
print("Maximum deviation: " + str(max_dev))
print("Deviation Acceptable: " + str(max_acceptable))

if max_dev > max_acceptable:
    print("Test failed")
else:
    print("Test passed")

# end the timer and display the duration
end = time.time()
print("Tests duration: " + str(end - start) + "s")
