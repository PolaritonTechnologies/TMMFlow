from FilterStack import FilterStack
import numpy as np
import pandas as pd
import time
import json


def execute_test():

    # Updated tests.py 05/09/2024
    json_path = "tests/TestStructure.json"

    # Read the JSON file
    with open(json_path, "r") as file:
        filter_json = json.load(file)

    # Create filter stack
    filter_stack = FilterStack(my_filter_dict=filter_json)

    # Wavelength is a numpy array of wavelengths in nm from 400 to 1200 nm spaced by 1 nm
    # Polar angles is a numpy array of polar angles in degrees from 0 to 89 degrees spaced by 1 degree
    # Azimuthal angles is a numpy array of azimuthal angles in degrees from 0 to 360 degrees spaced by 60 degree
    # Target type is a string that is either "r" or "t"
    # Polarization is a string that is either "s" or "p"

    wavelength = np.arange(400, 1200, 1)
    polar_angles = np.arange(0, 90, 1)
    # azimuthal_angles = np.arange(0, 360, 60)
    azimuthal_angles = np.arange(0, 1, 1)
    target_type = "r"
    polarization = "s"

    # Start a timer to evaluate how fast the calculation is executed using time library
    print("Starting test timer...")
    start = time.time()

    # Calculate the AR data
    filter_stack.calculate_ar_data(
        wavelength=wavelength,
        polar_angles=polar_angles,
        azimuthal_angles=azimuthal_angles,
        target_type=target_type,
        polarization=polarization,
        save_data=True,
        save_to_test=True,
    )

    # Stop the timer and output the duration inside the file test_performance.json where the first key is reference_test and the second is this_test
    end = time.time()
    duration = end - start
    print("Duration: ", duration, " seconds")

    # Before updating performance dictionary, we need to know that the test is passed in order to do so we compare test_to_compare.csv to this_test.csv
    # If this_test.csv does not exist then output that the reference test did not exist and return false
    # If this_test.csv does exist, compare the two files
    # If the two files are the same, output that the test passed and return true
    # If the two files are different, output that the test failed and return false

    # Check if this_test.csv exists
    try:
        with open("tests/this_test.csv", "r") as file:
            # load it with pandas
            this_test = pd.read_csv(file)
    except FileNotFoundError:
        print("Reference test does not exist, created with this execution")
        return False

    # Check if this_test.csv is the same as test_to_compare.csv
    with open("tests/test_reference.csv", "r") as file:
        test_to_compare = pd.read_csv(file)

    # Check if the two dataframes are the same with a certain threshold

    # Assuming this_test and test_to_compare are your DataFrames
    # Calculate the percentage of matching elements
    similarity_threshold = 0.995
    matching_elements = (this_test == test_to_compare).sum().sum()
    total_elements = this_test.size
    similarity_percentage = matching_elements / total_elements

    if similarity_percentage < similarity_threshold:
        return False

    # If the test passed, update the performance dictionary
    # If test_performance.json does not exist, create it and set the reference_test to the duration of this test
    try:
        with open("tests/test_performance.json", "r") as file:
            performance_dict = json.load(file)
    except FileNotFoundError:
        performance_dict = {"reference_test": duration}
    # If test_performance.json exists, set the this_test key to the duration of this test
    performance_dict["this_test"] = duration

    # Write the updated performance dictionary to test_performance.json
    with open("tests/test_performance.json", "w") as file:
        json.dump(performance_dict, file)
    return True


if __name__ == "__main__":
    test_result = execute_test()

    if test_result:
        print("Test passed")
    else:
        print("!!! Test failed !!!")
