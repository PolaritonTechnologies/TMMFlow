from FilterStack import FilterStack
import numpy as np
import pandas as pd
import time
import json


def execute_test(json_path, test_csv):

    with open(json_path, "r") as file:
        filter_json = json.load(file)

    # Create filter stack
    filter_stack = FilterStack(my_filter_dict=filter_json)

    # Wavelength is a numpy array of wavelengths in nm from 400 to 1200 nm spaced by 1 nm
    # Polar angles is a numpy array of polar angles in degrees from 0 to 89 degrees spaced by 1 degree
    # Azimuthal angles is a numpy array of azimuthal angles in degrees from 0 to 360 degrees spaced by 60 degree
    # Target type is a string that is either "r" or "t"
    # Polarization is a string that is either "s" or "p"

    wavelength = np.arange(
        filter_json["wavelengthMin"],
        filter_json["wavelengthMax"] + filter_json["wavelengthStep"],
        filter_json["wavelengthStep"],
    )
    polar_angles = np.arange(
        filter_json["polarAngleMin"],
        filter_json["polarAngleMax"] + filter_json["polarAngleStep"],
        filter_json["polarAngleStep"],
    )
    # azimuthal_angles = np.arange(0, 360, 60)
    azimuthal_angles = np.arange(
        filter_json["azimAngleMin"],
        filter_json["azimAngleMax"] + filter_json["azimAngleStep"],
        filter_json["azimAngleStep"],
    )
    target_type = filter_json["calculation_type"]
    polarization = filter_json["polarization"]

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
        save_to_test=False,
    )
    calculated_data = pd.DataFrame(np.array(filter_stack.stored_data)[0])

    # Stop the timer and print the duration of the simulation
    end = time.time()
    duration = end - start
    print("Duration: ", duration, " seconds")

    # Check if the reference exists and load it in from file
    try:
        with open(test_csv, "r") as file:
            # load it with pandas
            test_reference = pd.read_csv(file)
    except FileNotFoundError:
        print("Reference test does not exist: Please create one")

        return False

    # Now actually compare the just calculated data with the reference data
    similarity_threshold = 0.995
    matching_elements = (test_reference == calculated_data).sum().sum()
    total_elements = calculated_data.size
    similarity_percentage = matching_elements / total_elements

    print("Similarity percentage: ", similarity_percentage)

    if similarity_percentage < similarity_threshold:
        return False
    else:
        return True


if __name__ == "__main__":
    # Updated tests.py 05/09/2024
    json_paths = [
        "./src/tests/TestStructure.json",
        # "tests/best-39layers-AlPcCl-MgF2-860nm-885nm.json",
    ]

    # Compare the performance of the current test to the reference test (csv file)
    test_csvs = [
        "./src/tests/test_reference.csv",
        # "tests/best-39layers-AlPcCl-MgF2-860nm-885nm_reference.csv",
    ]

    # Read the JSON file
    for i in range(np.size(json_paths)):
        print("Running test ", i + 1, " of ", np.size(json_paths))

        test_result = execute_test(json_paths[i], test_csvs[i])

        if test_result:
            print("Test passed")
        else:
            print("!!! Test failed !!!")
