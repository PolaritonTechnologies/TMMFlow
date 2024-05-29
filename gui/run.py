from flask import (
    Flask,
    request,
    redirect,
    session,
    render_template,
    flash,
    jsonify,
    url_for,
    send_file,
)

# For asynchronous updates
from flask_socketio import SocketIO
from flask_cors import CORS

from flask.sessions import SecureCookieSessionInterface


# Plotting graphs
import plotly.graph_objects as go
import plotly.io as pio

# Threading computationally intensive tasks
import threading

# import subprocess
from werkzeug.utils import secure_filename

import os
import shutil
import time

import numpy as np
import pandas as pd
import json
import ast

from utility import allowed_file, generate_colors, get_available_materials, get_available_templates, is_number
from open_filter_converter import import_from_open_filter, export_to_open_filter
from FilterStack import FilterStack

class CustomSecureCookieSessionInterface(SecureCookieSessionInterface):
    def get_cookie_domain(self, app):
        return app.config.get('SESSION_COOKIE_DOMAIN')

    def get_cookie_path(self, app):
        return app.config.get('SESSION_COOKIE_PATH')

    def get_cookie_samesite(self, app):
        return app.config.get('SESSION_COOKIE_SAMESITE')

# Configure flask app
app = Flask(__name__)
app.session_interface = CustomSecureCookieSessionInterface()
app.config.update(
    SESSION_COOKIE_SAMESITE='Lax',  # or 'Strict' or 'None'
    SESSION_COOKIE_SECURE=True,
)

# Enable CORS for all routes which has to be done for the socketio to work,
# however, for the production server, it is important that cors_allowed_origins
# is only set to the server domain
CORS(app)
socketio = SocketIO(app, cors_allowed_origins="*")

# replace with your secret key
app.secret_key = "your secret key"
app.config["UPLOAD_FOLDER"] = "../src/temp/"
app.config["MATERIAL_FOLDER"] = "../materials/"
app.config["TEMPLATE_FOLDER"] = "../examples/"

# Global variables necessary for calculation (as data is shared)
my_filter = None
selected_file = None
default_file = "../examples/demo_test.json"

# Copy default file to temp folder and then use the copied file
shutil.copy(default_file, app.config["UPLOAD_FOLDER"] + default_file.split("/")[-1])
temp_default_file = app.config["UPLOAD_FOLDER"] + default_file.split("/")[-1]

# num_boxes = None
# colors = None
# heights = None
# unique_materials = None


##############################################
################## Routing ###################
##############################################


@app.route("/", methods=["GET", "POST"])
def stack():
    global my_filter
    global temp_default_file

    # If filter has not yet been selected, load the default filter and display
    # its structure, else
    if my_filter is None:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
            incoherent,
        ) = upload_file(temp_default_file)
    else:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
            incoherent
        ) = extract_filter_design(my_filter)

    # Specify the directory you want to search
    material_list = get_available_materials(app.config["MATERIAL_FOLDER"])
    template_list = get_available_templates(app.config["TEMPLATE_FOLDER"])

    default_values = {
        "num_boxes": num_boxes,
        "colors": colors,
        "heights": heights,
        "num_legend_items": num_legend_items,
        "unique_materials": unique_materials,
        "legend_colors": unique_colors,
        "file_label": selected_file,
        "available_materials": material_list,
        "available_templates": template_list,
    }

    # Add the entire filter definition to the default values to render it
    default_values.update(my_filter.filter_definition)
    temp_material_list = my_filter.structure_materials_by_user
    default_values["structure_materials"] = [
        (
            item.split("_")[0] + "_" + "_".join(item.split("_")[1:][::-1])
            if "_" in item
            else item
        )
        for item in temp_material_list
    ]
    default_values["structure_thicknesses"] = [
        np.flip(a).tolist() for a in my_filter.structure_thicknesses_by_user
    ]
    default_values["thickness_opt_allowed"] = my_filter.thickness_opt_allowed_by_user
    default_values["layer_switch_allowed"] = my_filter.layer_switch_allowed_by_user
    default_values["bounds"] = my_filter.bounds_by_user
    default_values["incoherent"] = my_filter.incoherent_by_user

    return render_template(
        "stack.html",
        default_values=default_values,
    )


@app.route("/simulate")
def simulate():
    global my_filter

    if np.all(my_filter.stored_data == None):
        dataPresent = False
    else:
        dataPresent = True

    default_values = {
        "mode": my_filter.filter_definition["calculation_type"],
        "startAngle": my_filter.filter_definition["polarAngleMin"],
        "endAngle": my_filter.filter_definition["polarAngleMax"],
        "stepAngle": my_filter.filter_definition["polarAngleStep"],
        "startWavelength": my_filter.filter_definition["wavelengthMin"],
        "endWavelength": my_filter.filter_definition["wavelengthMax"],
        "stepWavelength": my_filter.filter_definition["wavelengthStep"],
        "polarization": my_filter.filter_definition["polarization"],
        "azimuthalAngle": my_filter.filter_definition["azimAngleMin"],
        "dataPresent": dataPresent,
    }

    # Send the image URL back to the client
    return render_template(
        "simulate.html",
        default_values=default_values,
    )


@app.route("/optimize", methods=["GET", "POST"])
def optimize():
    global my_filter
    global temp_default_file

    if my_filter is None:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
            incoherent
        ) = upload_file(temp_default_file)
    else:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
            incoherent
        ) = extract_filter_design(my_filter)

    return render_template(
        "optimize.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=num_legend_items,
        unique_materials=unique_materials,
        legend_colors=unique_colors,
        incoherent= incoherent,
    )


@app.route("/materials")
def materials():
    # Get a list of all .csv files in the directory
    material_list = get_available_materials(app.config["MATERIAL_FOLDER"])

    return render_template("materials.html", available_materials=material_list)


##############################################
################## Stack ######################
##############################################


def identify_repeating_sequences(layers):
    layers = np.array(layers)
    sequences = [(layers[i], layers[i + 1]) for i in range(len(layers) - 1)]
    unique_sequences, counts = np.unique(sequences, return_counts=True, axis=0)
    repeating_sequences = unique_sequences[counts > 1]
    return repeating_sequences.tolist()


@app.route("/upload_file", methods=["POST"])
def upload_file(provided_filename=None):
    """
    Uploads a file to the server and processes it to generate filter stack
    representation.

    Args:
        defaultname (str, optional): The name of the default file to be used.
        Defaults to None.

    Returns:
        tuple or flask.Response: If defaultname is not None, returns a tuple
        containing the number of boxes, colors, heights, number of unique
        materials, unique materials, and unique colors.

        If defaultname is None, returns a flask.Response object with the
        rendered template.

    Raises:
        None
    """

    global selected_file
    global my_filter

    # Let the user upload a file to the server, if no file was chosen to be
    # uploaded, use the default file instead.
    if provided_filename == None:
        if "file" not in request.files:
            flash("No file part")
            return redirect(request.url)
        file = request.files["file"]
        if file.filename == "":
            flash("No selected file")
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            full_path = os.path.join(app.config["UPLOAD_FOLDER"], filename)
            file.save(full_path)

            if filename.split(".")[-1] == "json":
                # Standard .json format
                selected_file = full_path
            elif filename.split(".")[-1] == "ofp":
                # Open filter format (has to be translated to json first)
                json_path = import_from_open_filter(full_path)
                selected_file = json_path

    else:

        # If provided_filename is not none, use it as the selected file
        selected_file = provided_filename

    # Now check if the selected file is valid
    ###
    with open(selected_file, "r") as json_file_path:
        filter_definition_json = json.load(json_file_path)

    # Check that all materials in filter definition are available from the
    # materials DB
    temp_list_of_materials = [
        (
            a.split("_")[1:]
            if np.size(a.split("_")) > 1 and is_number(a.split("_")[0])
            else a
        )
        for a in np.unique(filter_definition_json["structure_materials"])
    ]
    list_of_materials_in_stack = np.unique(
        [
            item
            for sublist in temp_list_of_materials
            for item in (sublist if isinstance(sublist, list) else [sublist])
        ]
    )

    available_materials = get_available_materials(app.config["MATERIAL_FOLDER"])

    material_available = [
        material in available_materials for material in list_of_materials_in_stack
    ]

    # If the use selected replacements from the modal already
    if "replacements" in request.form:
        replacements = np.array(request.form["replacements"].split('"'))[1::4]

        # We need to reverse engineer, which of the following ones are not
        # available.
        # if filter_definition_json["substrate_material"] not in available_materials:
        # filter_definition_json["substrate_material"] = replacements[-1]
        # replacements = np.delete(replacements, -1)
        if filter_definition_json["incident_medium"] not in available_materials:
            filter_definition_json["incident_medium"] = replacements[-1]
            replacements = np.delete(replacements, -1)
        if filter_definition_json["exit_medium"] not in available_materials:
            filter_definition_json["exit_medium"] = replacements[-1]
            replacements = np.delete(replacements, -1)

        # Check that the size of replacements matches the size of materials that
        # were not found
        if np.size(replacements) == np.sum(np.invert(material_available)):
            # Iterate over all replacement elements, replace and save again to
            # file so that the correct file is loaded in
            for i in range(np.size(replacements)):
                temp_structure_materials = np.array(
                    filter_definition_json["structure_materials"]
                )
                temp_structure_materials[
                    np.array(filter_definition_json["structure_materials"])
                    == list_of_materials_in_stack[np.invert(material_available)][i]
                ] = replacements[i]
                filter_definition_json["structure_materials"] = (
                    temp_structure_materials.tolist()
                )

            # Now save again to file
            with open(selected_file, "w") as f:
                json.dump(filter_definition_json, f)

            material_available = True

    if (
        not np.all(material_available)
        or filter_definition_json["exit_medium"] not in available_materials
        or filter_definition_json["incident_medium"] not in available_materials
        # or filter_definition_json["substrate_material"] not in available_materials
    ):
        # This is a bit more complicated, as each of the materials could
        # potentially not be available
        unavailable_materials = []
        if not np.all(material_available):
            unavailable_materials = (
                unavailable_materials
                + list_of_materials_in_stack[np.invert(material_available)].tolist()
            )
        if filter_definition_json["exit_medium"] not in available_materials:
            unavailable_materials.append(filter_definition_json["exit_medium"])
        if filter_definition_json["incident_medium"] not in available_materials:
            unavailable_materials.append(filter_definition_json["incident_medium"])
        # if filter_definition_json["substrate_material"] not in available_materials:
            # unavailable_materials.append(filter_definition_json["substrate_material"])

        # open a modal to inform the user that some materials are not available
        # and break the loading (or better for the future: open a modal that
        # lets the user select the materials to choose instead)
        return jsonify(
            {
                "error": "Some materials are not available",
                "mismatched_materials": unavailable_materials,
                "available_materials": available_materials,
            }
        )

    ###

    # Create the filter and construct the filter stack representation
    my_filter = FilterStack(selected_file)
    (
        num_boxes,
        colors,
        heights,
        number_unique_materials,
        unique_materials,
        unique_colors,
        incoherent
    ) = extract_filter_design(my_filter)

    if provided_filename is not None:
        return (
            num_boxes,
            colors,
            heights,
            number_unique_materials,
            unique_materials,
            unique_colors,
            incoherent,
        )

    else:
        return jsonify({"redirect": url_for("stack")})
        # return render_template(
        # "stack.html",
        # num_boxes=num_boxes,
        # colors=colors,
        # heights=heights,
        # num_legend_items=number_unique_materials,
        # unique_materials=unique_materials,
        # legend_colors=unique_colors,
        # file_label=filename,
        # )




def extract_filter_design(filter):
    """
    Extracts the filter design information.

    Args:
        filter: The filter object containing the filter definition.

    Returns:
        A tuple containing the following information:
        - num_boxes: The number of boxes in the filter.
        - colors: An array of colors corresponding to each box.
        - heights: An array of heights corresponding to each box.
        - num_materials: The number of unique materials in the filter.
        - unique_materials: An array of unique materials in the filter.
        - unique_colors: An array of unique colors corresponding to each material.
    """

    # Now extract the number of layers to construct the number of boxes
    num_boxes = len(filter.filter_definition["structure_thicknesses"])
    unique_materials = np.unique(filter.filter_definition["structure_materials"])
    unique_colors = generate_colors(len(unique_materials))
    colors = np.empty(num_boxes, dtype=np.dtype("U7"))

    for i in range(len(unique_materials)):
        colors[
            np.array(filter.filter_definition["structure_materials"])
            == np.unique(filter.filter_definition["structure_materials"])[i]
        ] = unique_colors[i]
    
    ordered_thicknesses = np.array([
        filter.filter_definition["structure_thicknesses"][i] for i in filter.layer_order
    ])
    # Set the incoherent layer thicknesses to a specific value so that it does
    # not go out of control thick
    ordered_thicknesses[np.array(filter.filter_definition["incoherent"])] = 100

    heights = np.round(np.array(ordered_thicknesses) / sum(ordered_thicknesses) * 300)
    incoherent = filter.filter_definition["incoherent"]
    return (
        num_boxes,
        colors.tolist(),
        heights.tolist(),
        len(unique_materials),
        unique_materials.tolist(),
        unique_colors,
        incoherent
    )


@app.route("/save_json", methods=["POST"])
def save_json():
    global selected_file

    data = request.get_json()
    button = request.args.get("button")

    # Now the data has to be transformed into a json again so that it can be
    # saved to file. This is slightly annoying but I do not see a better of way
    # doing it right now.
    data_to_json = {}
    data_to_json["calculation_type"] = data[0].get("1").get("values")[0]
    data_to_json["polarization"] = data[1].get("1").get("values")[0]
    data_to_json["wavelengthMin"] = float(data[2].get("1").get("values")[0])
    data_to_json["wavelengthMax"] = float(data[2].get("2").get("values")[0])
    data_to_json["wavelengthStep"] = float(data[2].get("3").get("values")[0])
    data_to_json["polarAngleMin"] = float(data[3].get("1").get("values")[0])
    data_to_json["polarAngleMax"] = float(data[3].get("2").get("values")[0])
    data_to_json["polarAngleStep"] = float(data[3].get("3").get("values")[0])
    data_to_json["azimAngleMin"] = float(data[4].get("1").get("values")[0])
    data_to_json["azimAngleMax"] = float(data[4].get("2").get("values")[0])
    data_to_json["azimAngleStep"] = float(data[4].get("3").get("values")[0])
    # data_to_json["substrate_material"] = data[5].get("1").get("values")[0]
    # data_to_json["substrate_thickness"] = float(data[6].get("1").get("values")[0])
    data_to_json["incident_medium"] = data[5].get("1").get("values")[0]
    data_to_json["exit_medium"] = data[6].get("1").get("values")[0]

    layers = np.array(data[7].get("0").get("values"))
    number_of_columns = 7
    # Because of the way the structures are displayed (being different to the
    # order in the .json file), all entries concerning a layers must be flipped
    # here. Additionally, all entries containing muliple materials (e.g.
    # 5_SiO2_Ta2O5), have to be reverted, too (5_Ta2O5_SiO2) so that things make
    # sense. This is slighty confusing, however, and we should consider changing
    # the order in the .json file instead.
    temp_material_list = np.flip(layers[0::number_of_columns]).tolist()
    data_to_json["structure_materials"] = [
        (
            item.split("_")[0] + "_" + "_".join(item.split("_")[1:][::-1])
            if "_" in item
            else item
        )
        for item in temp_material_list
    ]
    temp_structure_thicknesses = np.flip(
        np.array(
            [ast.literal_eval(i) if "," in i else float(i) for i in layers[1::number_of_columns]],
            dtype=object,
        )
    ).tolist()
    data_to_json["structure_thicknesses"] = [
        np.flip(a).tolist() for a in temp_structure_thicknesses
    ]
    data_to_json["thickness_opt_allowed"] = np.flip(
        [s.lower() == "true" for s in layers[2::number_of_columns]]
    ).tolist()
    data_to_json["bounds"] = np.flip(
        [
            [float(x), float(y)] if y != "" else float(x)
            for x, y in zip(layers[3::number_of_columns].tolist(), layers[4::number_of_columns].tolist())
        ],
        axis=0,
    ).tolist()
    data_to_json["layer_switch_allowed"] = np.flip(
        [s.lower() == "true" for s in layers[5::number_of_columns]]
    ).tolist()
    data_to_json["incoherent"] = np.flip(
        [s.lower() == "true" for s in layers[6::number_of_columns]]
    ).tolist()

    targets = np.array(data[8].get("0").get("values"))

    no_columns_targets = 15
    data_to_json["targets_type"] = targets[0::no_columns_targets].tolist()
    data_to_json["targets_polarization"] = targets[1::no_columns_targets].tolist()
    data_to_json["targets_polar_angle"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(targets[2::no_columns_targets].tolist(), targets[3::no_columns_targets].tolist())
    ]
    data_to_json["polar_angle_steps"] = targets[4::no_columns_targets].astype(float).tolist()
    data_to_json["targets_azimuthal_angle"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(targets[5::no_columns_targets].tolist(), targets[6::no_columns_targets].tolist())
    ]
    data_to_json["azimuthal_angle_steps"] = targets[7::no_columns_targets].astype(float).tolist()
    data_to_json["targets_wavelengths"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(targets[8::no_columns_targets].tolist(), targets[9::no_columns_targets].tolist())
    ]
    data_to_json["wavelength_steps"] = targets[10::no_columns_targets].astype(float).tolist()
    data_to_json["targets_condition"] = targets[11::no_columns_targets].tolist()
    data_to_json["targets_value"] = targets[12::no_columns_targets].astype(float).tolist()
    data_to_json["targets_tolerance"] = targets[13::no_columns_targets].astype(float).tolist()
    data_to_json["targets_arithmetic"] = targets[14::no_columns_targets].astype(str).tolist()

    # If the file is not in the temp folder yet, do not allow for an overwrite
    if selected_file.split("/")[-2] != "temp":
        return jsonify({"error_message": "Cannot overwrite the default file."}), 200
    else:
        with open(selected_file, "w") as f:
            json.dump(data_to_json, f)

        # Reload filter for new data (the file path is still selected_file)
        upload_file(selected_file)
        return "", 302

        # if button == "download":
        #     # This is not working at the moment
        #     return send_file(selected_file, as_attachment=True)
        # else:


@app.route("/download")
def download_file():
    global my_filter

    file_ending = request.args.get("fileEnding")
    if file_ending == ".json":
        # Handle .json file ending
        path_to_file = selected_file
        return send_file(path_to_file, as_attachment=True)
    elif file_ending == ".ofp":
        # Handle .ofp file ending
        # First convert the cpp .json file to .ofp
        path_to_file = my_filter.json_file_path_cpp

        with open(path_to_file, "r") as input_file:
            input_dic = json.load(input_file)

        path_to_file = export_to_open_filter(
            input_dic, selected_file.split("/")[-1].split(".")[0] + "_converted"
        )

        return send_file(path_to_file, as_attachment=True)


@app.route("/download_current_optimum_file")
def download_current_optimum_file():
    global my_filter
    file_ending = request.args.get("fileEnding")
    path_to_json_file = os.path.join(
        app.config["UPLOAD_FOLDER"], "current_structure.json"
    )

    if file_ending == ".json":
        # Handle .json file ending
        return send_file(path_to_json_file, as_attachment=True)
    elif file_ending == ".ofp":
        # Handle .ofp file ending
        # First convert the .json file to .ofp
        input_file = path_to_json_file
        with open(input_file, "r") as input_file:
            input_dic = json.load(input_file)

        path_to_file = export_to_open_filter(input_dic, "current_structure_converted")

        return send_file(path_to_file, as_attachment=True)


@app.route("/reset_filter", methods=["POST"])
def reset_filter():
    global my_filter

    # Reset filter
    my_filter.reset()

    # Reload page
    return redirect("/")

@app.route('/start_new_design', methods=['POST'])
def start_new_design():
    filter_name = request.form.get('filter_name') + ".json"
    template = request.form.get('template') + ".json"

    new_file_path = os.path.join(app.config['UPLOAD_FOLDER'], filter_name)

    # Copy the template to the temp folder
    shutil.copy(os.path.join(app.config['TEMPLATE_FOLDER'], template), new_file_path)

    # Run the upload function with the new file
    upload_file(new_file_path)

    return '', 200


##############################################
############# Calculate and Plot #############
##############################################


@socketio.on("calculate_and_plot")
def calculate_and_plot(data):
    """
    Calculate AR data and plot
    """
    global my_filter
    wavelengths = np.arange(
        float(data["startWavelength"]),
        float(data["endWavelength"]) + float(data["stepWavelength"]),
        float(data["stepWavelength"]),
    )
    polar_angles = np.arange(
        float(data["startAngle"]),
        float(data["endAngle"]) + float(data["stepAngle"]),
        float(data["stepAngle"]),
    )
    # azimuthal_angles = np.arange(
    #     float(data["startAzimuthalAngle"]),
    #     float(data["endAzimuthalAngle"]) + float(data["stepAzimuthalAngle"]),
    #     float(data["stepAzimuthalAngle"]),
    # )
    target_type = data["mode"]
    polarization = "" if data["polarization"] == "None" else data["polarization"]

    my_filter.calculate_ar_data(
        wavelengths,
        polar_angles,
        azimuthal_angles=[float(data["azimuthalAngle"])],
        target_type=target_type,
        polarization=polarization,
    )
    calculated_data_df = my_filter.stored_data[0]

    # Create a Plotly figure using the calculated data
    angles = calculated_data_df.columns.to_numpy()
    wavelengths = calculated_data_df.index.to_numpy()
    color_values = calculated_data_df.to_numpy()

    # The layout is stored in simulate.html
    plotting_data = {
        "x": angles.tolist(),
        "y": wavelengths.tolist(),
        "z": color_values.tolist(),
    }

    # Emit the figure to the client
    socketio.emit("update_plot", plotting_data)

    """
    heatmap = go.Heatmap(
        x=angles,
        y=wavelengths,
        z=color_values,
        colorscale="Viridis",
    )

    fig = go.Figure(data=heatmap)

    # Convert the figure to JSON format
    fig_json = pio.to_json(fig)

    # Emit the figure in JSON format
    socketio.emit("update_plot", fig_json)
    """

    # calculated_data_df = pd.DataFrame()
    # for theta in polar_angles:
    #     calculated_data_df[theta] = my_filter.calculate_one_angle(
    #         float(data["startWavelength"]),
    #         float(data["endWavelength"]),
    #         float(data["stepWavelength"]),
    #         data["mode"],
    #         data["polarization"],
    #         theta,
    #         float(data["azimuthalAngle"]),
    #         True if data["generalCore"] == "on" else False,
    #     )

    #     # Create a Plotly figure using the calculated data
    #     angles = calculated_data_df.columns.to_numpy()
    #     wavelengths = calculated_data_df.index.to_numpy()
    #     color_values = calculated_data_df.to_numpy()

    #     # The layout is stored in simulate.html
    #     heatmap = go.Heatmap(
    #         x=angles,
    #         y=wavelengths,
    #         z=color_values,
    #         colorscale="Viridis",
    #     )

    #     fig = go.Figure(data=heatmap)

    #     # Convert the figure to JSON format
    #     fig_json = pio.to_json(fig)

    #     # Emit the figure in JSON format
    #     socketio.emit("update_plot", fig_json)

    # my_filter.stored_data = calculated_data_df


@socketio.on("plot")
def plot():
    """ """
    global my_filter

    calculated_data_df = my_filter.stored_data[0]

    # Create a Plotly figure using the calculated data
    angles = calculated_data_df.columns.to_numpy()
    wavelengths = calculated_data_df.index.to_numpy()
    color_values = calculated_data_df.to_numpy()

    # # The layout is stored in simulate.html
    # heatmap = go.Heatmap(
    #     x=angles,
    #     y=wavelengths,
    #     z=color_values,
    #     colorscale="Viridis",
    # )

    # fig = go.Figure(data=heatmap)

    # # Convert the figure to JSON format
    # fig_json = pio.to_json(fig)

    # # Emit the figure in JSON format
    # socketio.emit("update_plot", fig_json)

    # The layout is stored in simulate.html
    plotting_data = {
        "x": angles.tolist(),
        "y": wavelengths.tolist(),
        "z": color_values.tolist(),
    }

    # Emit the figure to the client
    socketio.emit("update_plot", plotting_data)


@socketio.on("plot_xy")
def handle_plot_xy(data):
    global my_filter
    x = my_filter.stored_data[0].index.to_list()

    # Get the y data corresponding to the x data
    y = my_filter.stored_data[0].loc[:, data["x"]].to_list()

    # Create a dictionary with the x and y data
    plot_data = {
        "x": x,
        "y": y,
        "name": str(data["x"])
        + "°, "
        + str(data["mode"])
        + ", "
        + str(data["polarization"]),
    }

    # Emit the data
    socketio.emit("update_xy_plot", plot_data)


@app.route("/download_data")
def download_data():
    global my_filter
    file_path = os.path.join(app.config["UPLOAD_FOLDER"], "simulated_data.csv")

    x = my_filter.stored_data[0]
    x.to_csv(file_path, sep="\t")

    return send_file(file_path, as_attachment=True)


##############################################
################ Optimization  ###############
##############################################


@socketio.on("start_optimization")
def start_optimization(data):
    """ """
    global my_filter

    # my_filter.perform_optimisation(
    # data["optimizationMethod"],
    # save_optimized_to_file=True,
    # )
    thread = threading.Thread(
        target=my_filter.perform_optimisation, args=(data["optimizationMethod"],)
    )
    thread.start()

    # Some helper variables for plotting
    i = 0
    current_optimization_method = data["optimizationMethod"][i]
    colors = generate_colors(len(data["optimizationMethod"]))
    traces = []
    traces.append(
        {"x": [], "y": [], "color": colors[i], "name": current_optimization_method}
    )

    time.sleep(0.1)

    while not my_filter.stop_flag:
        if my_filter.optimization_method != current_optimization_method:
            # Plot an additional line in the merit graph indicating the switch
            # to a new method
            i += 1
            current_optimization_method = data["optimizationMethod"][i]

            # Start a new trace
            traces.append(
                {
                    "x": [],
                    "y": [],
                    "color": colors[i],
                    "name": current_optimization_method,
                }
            )

        # Add the current data to the current trace
        traces[-1]["x"].append(my_filter.optimum_iteration)
        traces[-1]["y"].append(my_filter.optimum_merit)

        # Emit the traces to the client
        socketio.emit("update_merit_graph", traces)

        # Update the stack representation
        (
            num_boxes,
            colors,
            heights,
            number_unique_materials,
            unique_materials,
            unique_colors,
            incoherent
        ) = extract_filter_design(my_filter)

        # Package the values into a dictionary
        filter_representation = {
            "num_boxes": num_boxes,
            "colors": colors,
            "heights": heights,
            "number_unique_materials": number_unique_materials,
            "unique_materials": unique_materials,
            "unique_colors": unique_colors,
            "incoherent": incoherent,
        }

        # Convert the dictionary to a JSON string
        filter_json = json.dumps(filter_representation)

        # Emit the JSON string
        socketio.emit(
            "update_filter_representation",
            {
                "filter_json": filter_json,
                "iterations": my_filter.iteration_no,
                "merit": np.round(my_filter.last_merit, 1),
            },
        )

        time.sleep(0.5)


@socketio.on("stop_optimization")
def stop_optimization():
    """ """
    global my_filter

    my_filter.stop_flag = True


##############################################
################# Materials  #################
##############################################


@socketio.on("get_material_data")
def get_material_data(data):
    material = data["material"]
    # Load file and get data
    df = pd.read_csv("../materials/" + material + ".csv", skiprows=1, sep="\t")

    # Convert the first column to a list and all following columns to a list of lists
    data = {
        "x": df.iloc[:, 0].tolist(),
        "y": df.iloc[:, 1:].values.T.tolist(),
        "name": df.columns[1:].tolist(),
    }

    socketio.emit("material_data", data)


@app.route("/download_material")
def download_material():
    file_name = request.args.get("fileName")
    path_to_file = app.config["MATERIAL_FOLDER"] + file_name + ".csv"

    return send_file(path_to_file, as_attachment=True)


@app.route("/upload_material", methods=["POST"])
def upload_material():
    if "file" not in request.files:
        return "No file part", 400
    file = request.files["file"]
    if file.filename == "":
        return "No selected file", 400
    if file:
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config["UPLOAD_FOLDER"], filename))

        # Load file and get data
        df = pd.read_csv(
            os.path.join(app.config["UPLOAD_FOLDER"], filename), skiprows=1, sep="\t"
        )

        # Convert the first column to a list and all following columns to a list of lists
        data = {
            "x": df.iloc[:, 0].tolist(),
            "y": df.iloc[:, 1:].values.T.tolist(),
            "name": df.columns[1:].tolist(),
        }

        # Emit the data
        socketio.emit("material_data", data)

        return filename, 200


@app.route("/accept", methods=["POST"])
def accept_file():
    filename = request.form.get("filename")
    shutil.move(
        os.path.join(app.config["UPLOAD_FOLDER"], filename),
        os.path.join(app.config["MATERIAL_FOLDER"], filename),
    )
    return "File accepted", 200


@app.route("/reject", methods=["POST"])
def reject_file():
    filename = request.form.get("filename")
    os.remove(os.path.join(app.config["UPLOAD_FOLDER"], filename))
    return "File rejected", 200


if __name__ == "__main__":
    socketio.run(app)
