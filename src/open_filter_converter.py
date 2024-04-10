import json
import numpy as np


def import_from_open_filter(input_text):

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
        "targets_wavelengths": [],
        "wavelength_steps": [],
        "targets_tolerance": [],
        "bounds": [],
        "add_layers": False,
        "nb_added_layers": 0,
        "added_materials": [],
        "added_layer_bounds": [],
    }

    # Go through the open filter file to populate the dictionary
    for idx, line in enumerate(lines):
        if line.startswith("FrontLayer:"):
            parts = line.split(" ")
            if parts[1] == "Ta2O5_Koeln_AM":
                data["structure_materials"].append("Ta2O5")
            elif parts[1] == "SiO2_Koeln_AM":
                data["structure_materials"].append("SiO2")
            else:
                data["structure_materials"].append(parts[1])
            data["structure_thicknesses"].append(float(parts[2]))
        if line.startswith("RefineThickness:"):
            parts = line.split(" ")
            if parts[1] == "1":
                data["thickness_opt_allowed"].append(True)
                # To be changed later according to what is needed
                data["layer_switch_allowed"].append(False)
            elif parts[1] == "0":
                data["thickness_opt_allowed"].append(False)
                # To be changed later according to what is needed
                data["layer_switch_allowed"].append(False)
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
    # Fill other fields by default
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
    if data["wavelengthMin"] is None:
        data["wavelengthMin"] = 400.0
    if data["wavelengthMax"] is None:
        data["wavelengthMax"] = 800.0
    if data["wavelengthStep"] is None:
        data["wavelengthStep"] = 1.0
    if data["polar_angle_steps"] == []:
        data["polar_angle_steps"] = list(np.ones_like(data["targets_tolerance"]))
    if data["targets_azimuthal_angle"] == []:
        data["targets_azimuthal_angle"] = list(np.zeros_like(data["targets_tolerance"]))
    if data["azimuthal_angle_steps"] == []:
        data["azimuthal_angle_steps"] = list(np.ones_like(data["targets_tolerance"]))
    if data["bounds"] == []:
        data["bounds"] = [[0, 200]] * len(data["structure_materials"])
    print(data)
    # Write to JSON file
    with open("converted_open_filter.json", "w") as output_file:
        json.dump(data, output_file)


def export_to_open_filter(dictionary_input):

    to_export = ""

    header = """Version: 1.1.1
Comment:
End
Filter:
    Substrate: substrate_24mm_koeln 1000000.000000
    FrontMedium: void
    BackMedium: void
    CenterWavelength: 450.000000
    WavelengthRange: 300.000000 1000.000000 1.000000
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

    for i, el in enumerate(dictionary_input["structure_materials"]):

        if el == "Ta2O5":

            material_block = f"""   FrontLayer: Ta2O5_Koeln_AM {dictionary_input['structure_thicknesses'][i]}
    RefineThickness: {int(dictionary_input['thickness_opt_allowed'][i])}
"""
            to_export = to_export + material_block

        elif el == "SiO2":

            material_block = f"""   FrontLayer: SiO2_Koeln_AM {dictionary_input['structure_thicknesses'][i]}
    RefineThickness: {int(dictionary_input['thickness_opt_allowed'][i])}
"""
            to_export = to_export + material_block

        else:

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
    Angle: {dictionary_input['targets_polar_angle'][i]}
    Polarization: {dictionary_polarization[dictionary_input['targets_polarization'][i]]}
    Direction: Normal
    From: {dictionary_input['targets_wavelengths'][i][0]}
    To: {dictionary_input['targets_wavelengths'][i][1]}
    By: {dictionary_input['wavelength_steps'][i]}
    Points: 
        {dictionary_input['targets_wavelengths'][i][0]}\t{dictionary_input['targets_value'][i]}\t{dictionary_input['targets_tolerance'][i]}
        {dictionary_input['targets_wavelengths'][i][1]}\t{dictionary_input['targets_value'][i]}\t{dictionary_input['targets_tolerance'][i]}
    End"""
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
End"""

            if dictionary_input["targets_condition"][i] != "":
                target_block = (
                    target_block
                    + f"""Inequality: {dictionary_condition[dictionary_input['targets_condition'][i]]}
    End
    """
                )
        to_export = to_export + target_block

    # Write to text file
    with open("converted_open_filter.ofp", "w") as output_file:
        output_file.write(to_export)


if True:

    input_file = "test_new_core.ofp"
    with open(input_file, "r") as input_file:
        input_text = input_file.read()
    print(import_from_open_filter(input_text))

if False:

    input_file = "current_structure.json"
    with open(input_file, "r") as input_file:
        input_dic = json.load(input_file)
        print(input_dic)
    print(export_to_open_filter(input_dic))
