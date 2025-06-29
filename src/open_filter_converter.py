import json
import numpy as np
import os
import re


def import_from_open_filter(input_file):
    with open(input_file, "r") as input_file:
        input_text = input_file.read()

    # Initial read of the import open filter project
    lines = input_text.split("\n")

    # Empty dictionary gets instantiated here
    data = {
        "calculation_type": "",
        "polarAngleMin": None,
        "polarAngleMax": None,
        "polarAngleStep": None,
        "azimAngleMin": None,
        "azimAngleMax": None,
        "azimAngleStep": None,
        "wavelengthMin": None,
        "wavelengthMax": None,
        "wavelengthStep": None,
        "structure_materials": [],
        "structure_thicknesses": [],
        "polarization": 0.5,
        "thickness_opt_allowed": [],
        "layer_switch_allowed": [],
        "targets_type": [],
        "targets_polarization": [],
        "targets_condition": [],
        "targets_value": [],
        "targets_polar_angle": [],
        "polar_angle_steps": [],
        "targets_azimuthal_angle": [],
        "azimuthal_angle_steps": [],
        "incoherent": [],
        "targets_wavelengths": [],
        "wavelength_steps": [],
        "targets_tolerance": [],
        "targets_arithmetic": [],
        "bounds": [],
        "incident_medium": None,
        "exit_medium": None,
        "core_selection": None,
    }

    # Go through the open filter file to populate the dictionary
    front_layer = True
    for idx, line in enumerate(lines):
        # Remove leading or trailing white space
        line = line.strip()
        parts = re.split(r"\s+", line)  # Split by any whitespace

        #### General ####
        if line.startswith("FrontMedium:"):
            if parts[1] == "void":
                data["incident_medium"] = "Air"
            else:
                data["incident_medium"] = str(parts[1])
        elif line.startswith("BackMedium:"):
            if parts[1] == "void":
                data["exit_medium"] = "Air"
            else:
                data["exit_medium"] = str(parts[1])
        elif line.startswith("WavelengthRange:"):
            data["wavelengthMin"] = float(parts[1])
            data["wavelengthMax"] = float(parts[2])
            data["wavelengthStep"] = float(parts[3])

        #### Materials ####
        elif line.startswith("Substrate:"):
            data["structure_materials"].append(parts[1])
            data["structure_thicknesses"].append(float(parts[2]))
            data["thickness_opt_allowed"].append(False)
            data["layer_switch_allowed"].append(False)
            data["incoherent"].append(True)
        elif line.startswith("FrontLayer:"):
            front_layer = True
            data["structure_materials"].append(parts[1])
            data["structure_thicknesses"].append(float(parts[2]))
            data["incoherent"].append(False)
        elif line.startswith("BackLayer:"):
            front_layer = False
            data["structure_materials"].insert(0, parts[1])
            data["structure_thicknesses"].insert(0, float(parts[2]))
            data["incoherent"].insert(0, False)
        elif line.startswith("RefineThickness:"):
            if parts[1] == "1":
                if front_layer:
                    data["thickness_opt_allowed"].append(True)
                    data["layer_switch_allowed"].append(False)
                else:
                    data["thickness_opt_allowed"].insert(0, True)
                    data["layer_switch_allowed"].insert(0, False)
            elif parts[1] == "0":
                if front_layer:
                    data["thickness_opt_allowed"].append(False)
                    data["layer_switch_allowed"].append(False)
                else:
                    data["thickness_opt_allowed"].insert(0, False)
                    data["layer_switch_allowed"].insert(0, False)

        #### Targets ####
        elif line.startswith("Kind:"):
            if parts[1] == "TransmissionSpectrum":
                data["targets_type"].append("t")
            elif parts[1] == "AbsorptionSpectrum":
                data["targets_type"].append("a")
            elif parts[1] == "ReflectionSpectrum":
                data["targets_type"].append("r")
            elif parts[1] == "Transmission":
                data["targets_type"].append("t")
            elif parts[1] == "Absorption":
                data["targets_type"].append("a")
            elif parts[1] == "Reflection":
                data["targets_type"].append("r")
        elif line.startswith("Angle:"):
            data["targets_polar_angle"].append(float(parts[1]))
        elif line.startswith("Polarization:"):
            if float(parts[1]) == 45.0:
                data["targets_polarization"].append(0.5)
            elif float(parts[1]) == 90.0:
                data["targets_polarization"].append(1.0)
            elif float(parts[1]) == 0.0:
                data["targets_polarization"].append(0.0)
        elif line.startswith("From:"):
            temp_array = [float(parts[1])]
            parts_next_line = re.split(r"\s+", lines[idx + 1].strip())
            temp_array.append(float(parts_next_line[-1]))
            data["targets_wavelengths"].append(temp_array)
        elif line.startswith("Wavelength:"):
            data["targets_wavelengths"].append(float(parts[1]))
        elif line.startswith("By:"):
            data["wavelength_steps"].append(float(parts[1]))
        elif line.startswith("Points:"):
            parts_next_line = re.split(r"\s+", lines[idx + 1].strip())
            data["targets_value"].append(float(parts_next_line[1]))
            data["targets_tolerance"].append(float(parts_next_line[2]))
            data["targets_condition"].append("=")
        elif line.startswith("Inequality:"):
            if parts[1] == "larger":
                data["targets_condition"].append(">")
            elif parts[1] == "smaller":
                data["targets_condition"].append("<")
        elif line.startswith("Value:"):
            data["targets_value"].append(float(parts[1]))
        elif line.startswith("Tolerance:"):
            data["targets_tolerance"].append(float(parts[1]))

    #### Populate Fields that do not exist in OpenFilters ####
    if data["calculation_type"] == "":
        data["calculation_type"] = "t"
    if data["polarAngleMin"] is None:
        data["polarAngleMin"] = 0.0
    if data["polarAngleMax"] is None:
        data["polarAngleMax"] = 80.0
    if data["polarAngleStep"] is None:
        data["polarAngleStep"] = 1.0
    if data["azimAngleMin"] is None:
        data["azimAngleMin"] = 0.0
    if data["azimAngleMax"] is None:
        data["azimAngleMax"] = 0.0
    if data["azimAngleStep"] is None:
        data["azimAngleStep"] = 1.0
    if data["polar_angle_steps"] == []:
        data["polar_angle_steps"] = list(np.ones_like(data["targets_tolerance"]))
    if data["targets_azimuthal_angle"] == []:
        data["targets_azimuthal_angle"] = list(np.zeros_like(data["targets_tolerance"]))
    if data["azimuthal_angle_steps"] == []:
        data["azimuthal_angle_steps"] = list(np.ones_like(data["targets_tolerance"]))
    if data["targets_arithmetic"] == []:
        data["targets_arithmetic"] = list(
            np.arange(1, len(data["targets_tolerance"]) + 1, dtype=int).astype(str)
        )
    if data["bounds"] == []:
        data["bounds"] = [[0, 500]] * len(data["structure_materials"])
    if data["core_selection"] is None:
        data["core_selection"] = "fast"

    # Write to JSON file
    saving_path = input_file.name[:-3] + "json"
    with open(saving_path, "w") as output_file:
        json.dump(data, output_file)

    return saving_path


def export_to_open_filter(dictionary_input, file_name):

    if np.sum(dictionary_input["incoherent"]) > 1:
        message = (
            "Export to open filters is only possible with exactly one incoherent layer."
        )
        print(message)
        return message

    # Currently no support for exporting arbitrary polarizations to open filters (although it is possible)
    if np.any(
        [
            target_polarization not in [0.0, 1.0, 0.5]
            for target_polarization in dictionary_input["targets_polarization"]
        ]
    ):
        message = "Only polarizations of 0.0, 0.5, and 1.0 are supported for open filters export at this point"
        print(message)
        return message

    to_export = ""
    substrate_index = dictionary_input["incoherent"].index(True)
    substrate = dictionary_input["structure_materials"][(substrate_index)]
    substrate_thickness = dictionary_input["structure_thicknesses"][substrate_index]

    header = f"""Version: 1.1.1
    Comment:
    End
    Filter:
		Substrate: {substrate} {substrate_thickness}
		FrontMedium: {dictionary_input['incident_medium'] if dictionary_input == 'Air' else 'void'}
		BackMedium: {dictionary_input['exit_medium'] if dictionary_input == 'Air' else 'void'}
		CenterWavelength: 450.000000
		WavelengthRange: {dictionary_input['wavelengthMin']} {dictionary_input['wavelengthMax']} {dictionary_input['wavelengthStep']}
		DontConsiderSubstrate: 0
		StepSpacing: 0.010000
		MinimumThickness: 0.000000
		Illuminant: CIE-D65
		Observer: CIE-1931
		ConsiderBackside: 1
		EllipsometerType: 1
		DeltaMin: -90.000000
		ConsiderBacksideOnMonitoring: 1
		MonitoringEllipsometerType: 1
		MonitoringDeltaMin: -90.000000
		MonitoringSublayerThickness: 1.000000
    """

    to_export = to_export + header

    # Front layers
    for i, material in enumerate(
        dictionary_input["structure_materials"][substrate_index + 1 :]
    ):
        current_index = substrate_index + 1 + i

        material_block = f"""	FrontLayer: {material} {dictionary_input['structure_thicknesses'][current_index]}
	RefineThickness: {int(dictionary_input['thickness_opt_allowed'][current_index])}
"""
        to_export = to_export + material_block

    # Back layers

    if substrate_index != 0:
        back_materials = dictionary_input["structure_materials"][:substrate_index]
        back_materials.reverse()
        for i, material in enumerate(back_materials):
            current_index = len(back_materials) - 1 - i

            material_block = f"""	BackLayer: {material} {dictionary_input['structure_thicknesses'][current_index]}
	RefineThickness: {int(dictionary_input['thickness_opt_allowed'][current_index])}
    """
            to_export = to_export + material_block

    to_export = to_export + "End\n"

    # Define the dictionaries
    dictionary_polarization = {1.0: 0.0, 0.0: 90.0, 0.5: 45.0}

    dictionary_type = {
        "t": "TransmissionSpectrum",
        "r": "ReflectionSpectrum",
        "a": "AbsorptionSpectrum",
    }

    dictionary_type_single = {
        "t": "Transmission",
        "r": "Reflection",
        "a": "Absorption",
    }

    dictionary_condition = {
        "=": "",
        ">": "Inequality: larger",
        "<": "Inequality: smaller",
    }

    # Helper functions
    def generate_target_block(dictionary_input, i, angle=None):
        target_block = f"""Target:
        Kind: {dictionary_type[dictionary_input['targets_type'][i]]}
        Angle: {angle if angle is not None else dictionary_input['targets_polar_angle'][i]}
        Polarization: {dictionary_polarization[dictionary_input['targets_polarization'][i]]}
        Direction: Normal
        From: {dictionary_input['targets_wavelengths'][i][0]}
        To: {dictionary_input['targets_wavelengths'][i][1]}
        By: {dictionary_input['wavelength_steps'][i]}
        Points: 
            {dictionary_input['targets_wavelengths'][i][0]}        {dictionary_input['targets_value'][i]}        {dictionary_input['targets_tolerance'][i]}
            {dictionary_input['targets_wavelengths'][i][1]}        {dictionary_input['targets_value'][i]}        {dictionary_input['targets_tolerance'][i]}
        End
    """
        if dictionary_input["targets_condition"][i] != "":
            target_block += f"""
        {dictionary_condition[dictionary_input['targets_condition'][i]]}
    End
    """
        return target_block

    def generate_single_target_block(dictionary_input, i, angle=None):
        target_block = f"""Target:
        Kind: {dictionary_type_single[dictionary_input['targets_type'][i]]}
        Angle: {angle if angle is not None else dictionary_input['targets_polar_angle'][i]}
        Polarization: {dictionary_polarization[dictionary_input['targets_polarization'][i]]}
        Direction: Normal
        Wavelength: {dictionary_input['targets_wavelengths'][i]}
        Value: {dictionary_input['targets_value'][i]}
        Tolerance: {dictionary_input['targets_tolerance'][i]}
    """
        if dictionary_input["targets_condition"][i] != "":
            target_block += f"""
        {dictionary_condition[dictionary_input['targets_condition'][i]]}
    End
    """
        return target_block

    # Main loop to generate the export string
    for i, el in enumerate(dictionary_input["targets_type"]):
        target_block = ""

        # Handle polar angles
        if (
            isinstance(dictionary_input["targets_polar_angle"][i], list)
            and len(dictionary_input["targets_polar_angle"][i]) == 2
        ):
            polar_angle_min = dictionary_input["targets_polar_angle"][i][0]
            polar_angle_max = dictionary_input["targets_polar_angle"][i][1]
            polar_angle_step = dictionary_input["polar_angle_steps"][i]

            polar_angles = np.arange(
                polar_angle_min,
                polar_angle_max + polar_angle_step,
                polar_angle_step,
            )

            for angle in polar_angles:
                if isinstance(dictionary_input["targets_wavelengths"][i], list):
                    target_block += generate_target_block(dictionary_input, i, angle)
                else:
                    target_block += generate_single_target_block(
                        dictionary_input, i, angle
                    )
        else:
            if isinstance(dictionary_input["targets_wavelengths"][i], list):
                target_block = generate_target_block(dictionary_input, i)
            else:
                target_block = generate_single_target_block(dictionary_input, i)

        to_export += target_block

    # Write to text file
    path_to_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "./temp/" + file_name + ".ofp"
    )

    with open(
        path_to_file,
        "w",
    ) as output_file:
        output_file.write(to_export)

    return path_to_file


def convert_csv_from_open_filter(file_path):
    import pandas as pd

    temp = pd.read_csv(file_path, skiprows=1)
    wavelength = temp["wavelength"].unique().reshape(1, -1)
    absorption = temp["absorption"].to_numpy()

    absorption_reshaped = np.reshape(absorption, [7, int(len(absorption) / 7)])
    new_columns = np.array(["wavelength", "0", "15", "30", "45", "60", "75", "89"])

    temp_save_df = pd.DataFrame(
        np.concatenate((wavelength, absorption_reshaped)).T, columns=new_columns
    )
    temp_save_df.to_csv(file_path + "_converted.csv", index=False, sep="\t")


# if True:

#     input_file = "test_new_core.ofp"
#     with open(input_file, "r") as input_file:
#         input_text = input_file.read()
#     print(import_from_open_filter(input_text))

# if True:

#     input_file = "../examples/current_structure_bestAlPCCl.json"
#     with open(input_file, "r") as input_file:
#         input_dic = json.load(input_file)
#         print(input_dic)
#     print(export_to_open_filter(input_dic))
