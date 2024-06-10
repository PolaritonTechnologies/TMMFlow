import json
import numpy as np
import os


def import_from_open_filter(input_file):
    with open(input_file, "r") as input_file:
        input_text = input_file.read()

    # Initial read of the import open filter project
    lines = [line.replace("\t", "") for line in input_text.split("\n")]

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
        "polarization": "",
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
        #### General ####
        if line.startswith("FrontMedium:"):
            parts = line.split(" ")
            if parts[1] == "void":
                data["incident_medium"] = "Air"
            else:
                data["incident_medium"] = str(parts[1])
        if line.startswith("BackMedium:"):
            parts = line.split(" ")
            if parts[1] == "void":
                data["exit_medium"] = "Air"
            else:
                data["exit_medium"] = str(parts[1])
        if line.startswith("WavelengthRange:"):
            parts = line.split(" ")
            data["wavelengthMin"] = float(parts[1])
            data["wavelengthMax"] = float(parts[2])
            data["wavelengthStep"] = float(parts[3])

        #### Materials ####
        if line.startswith("Substrate:"):
            parts = line.split(" ")
            data["structure_materials"].append(parts[1])
            data["structure_thicknesses"].append(float(parts[2]))
            data["thickness_opt_allowed"].append(False)
            data["layer_switch_allowed"].append(False)
            data["incoherent"].append(True)
        if line.startswith("FrontLayer:"):
            front_layer = True
            parts = line.split(" ")
            data["structure_materials"].append(parts[1])
            data["structure_thicknesses"].append(float(parts[2]))
            data["incoherent"].append(False)
        if line.startswith("BackLayer:"):
            front_layer = False
            parts = line.split(" ")
            data["structure_materials"].insert(0, parts[1])
            data["structure_thicknesses"].insert(0, float(parts[2]))
            data["incoherent"].insert(0, False)
        if line.startswith("RefineThickness:"):
            parts = line.split(" ")
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
        if line.startswith("Kind:"):
            parts = line.split(" ")
            if parts[1] == "TransmissionSpectrum":
                data["targets_type"].append("t")
            if parts[1] == "AbsorptionSpectrum":
                data["targets_type"].append("a")
            if parts[1] == "ReflectionSpectrum":
                data["targets_type"].append("r")
            if parts[1] == "Transmission":
                data["targets_type"].append("t")
            if parts[1] == "Absorption":
                data["targets_type"].append("a")
            if parts[1] == "Reflection":
                data["targets_type"].append("r")
        if line.startswith("Angle:"):
            parts = line.split(" ")
            data["targets_polar_angle"].append(float(parts[1]))
        if line.startswith("Polarization:"):
            parts = line.split(" ")
            print(parts)
            if float(parts[1]) == 45.0:
                data["targets_polarization"].append("")
            if float(parts[1]) == 90.0:
                data["targets_polarization"].append("s")
            if float(parts[1]) == 0.0:
                data["targets_polarization"].append("p")
        if line.startswith("From:"):
            parts = line.split(" ")
            temp_array = []
            temp_array.append(float(parts[1]))
            parts_next_line = lines[idx + 1].split(" ")
            temp_array.append(float(parts_next_line[1]))
            data["targets_wavelengths"].append(temp_array)
        if line.startswith("Wavelength:"):
            parts = line.split(" ")
            data["targets_wavelengths"].append(float(parts[1]))
        if line.startswith("By:"):
            parts = line.split(" ")
            data["wavelength_steps"].append(float(parts[1]))
        if line.startswith("Points:"):
            parts_next_line = lines[idx + 1].split(" ")
            parts_next_line = [item for item in parts_next_line if item != ""]
            data["targets_value"].append(float(parts_next_line[1]))
            data["targets_tolerance"].append(float(parts_next_line[-1]))
            data["targets_condition"].append("=")
        if line.startswith("Inequality:"):
            parts = line.split(" ")
            if parts[1] == "larger":
                data["targets_condition"].append(">")
            if parts[1] == "smaller":
                data["targets_condition"].append("<")
        if line.startswith("Value:"):
            parts = line.split(" ")
            data["targets_value"].append(float(parts[1]))
        if line.startswith("Tolerance:"):
            parts = line.split(" ")
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

    to_export = ""

    header = f"""Version: 1.1.1
    Comment:
    End
    Filter:
        Substrate: {dictionary_input['structure_materials'][0]} {dictionary_input['structure_thicknesses'][0]}
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

    for i, el in enumerate(dictionary_input["structure_materials"][1:]):

        material_block = f"""   FrontLayer: {dictionary_input['structure_materials'][i]} {dictionary_input['structure_thicknesses'][i]}
    RefineThickness: {int(dictionary_input['thickness_opt_allowed'][i])}
"""
        to_export = to_export + material_block

    to_export = to_export + "End\n"

    dictionary_polarization = {"s": 0.0, "p": 90.0, "": 45.0}
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

    for i, el in enumerate(dictionary_input["targets_type"]):

        if isinstance(dictionary_input["targets_wavelengths"][i], list):
            target_block = f"""Target:
    Kind: {dictionary_type[dictionary_input['targets_type'][i]]}
    Angle: {dictionary_input['targets_polar_angle'][i] if np.size(dictionary_input['targets_polar_angle'][i]) == 1 else dictionary_input['targets_polar_angle'][i][0]}
    Polarization: {dictionary_polarization[dictionary_input['targets_polarization'][i]]}
    Direction: Normal
    From: {dictionary_input['targets_wavelengths'][i][0]}
    To: {dictionary_input['targets_wavelengths'][i][1]}
    By: {dictionary_input['wavelength_steps'][i]}
    Points: 
        {dictionary_input['targets_wavelengths'][i][0]}\t{dictionary_input['targets_value'][i]}\t{dictionary_input['targets_tolerance'][i]}
        {dictionary_input['targets_wavelengths'][i][1]}\t{dictionary_input['targets_value'][i]}\t{dictionary_input['targets_tolerance'][i]}
    End
"""
            if dictionary_input["targets_condition"][i] != "":
                target_block = (
                    target_block
                    + f"""
        {dictionary_condition[dictionary_input['targets_condition'][i]]}
End
"""
                )
        else:
            target_block = f"""Target:
    Kind: {dictionary_type_single[dictionary_input['targets_type'][i]]}
    Angle: {dictionary_input['targets_polar_angle'][i]}
    Polarization: {dictionary_polarization[dictionary_input['targets_polarization'][i]]}
    Direction: Normal
    Wavelength: {dictionary_input['targets_wavelengths'][i]}
    Value: {dictionary_input['targets_value'][i]}
    Tolerance: {dictionary_input['targets_tolerance'][i]}
    """

            if dictionary_input["targets_condition"][i] != "":
                target_block = (
                    target_block
                    + f"""{dictionary_condition[dictionary_input['targets_condition'][i]]}
End
"""
                )
        to_export = to_export + target_block

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
