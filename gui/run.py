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

from utility import (
    allowed_file,
    generate_colors,
)
from open_filter_converter import import_from_open_filter, export_to_open_filter
from FilterStack import FilterStack

# Configure flask app
app = Flask(__name__)

# Enable CORS for all routes which has to be done for the socketio to work,
# however, for the production server, it is important that cors_allowed_origins
# is only set to the server domain
CORS(app)
socketio = SocketIO(app, cors_allowed_origins="*")

# replace with your secret key
app.secret_key = "your secret key"
app.config["UPLOAD_FOLDER"] = "../src/temp/"
app.config["MATERIAL_FOLDER"] = "../materials/"

# Global variables necessary for calculation (as data is shared)
my_filter = None
selected_file = None
default_file = "../examples/demo_test.json"

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
    global default_file

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
        ) = upload_file(default_file)
    else:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
        ) = extract_filter_design(my_filter)

    # Specify the directory you want to search
    directory = "../materials/"

    material_list = [
        os.path.splitext(f)[0]
        for f in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, f)) and f.endswith(".csv")
    ]

    default_values = {
        "num_boxes": num_boxes,
        "colors": colors,
        "heights": heights,
        "num_legend_items": num_legend_items,
        "unique_materials": unique_materials,
        "legend_colors": unique_colors,
        "file_label": selected_file,
        "available_materials": material_list,
    }

    # Add the entire filter definition to the default values to render it
    default_values.update(my_filter.filter_definition)
    default_values["structure_materials"] = my_filter.structure_materials_by_user
    default_values["structure_thicknesses"] = my_filter.structure_thicknesses_by_user
    default_values["thickness_opt_allowed"] = my_filter.thickness_opt_allowed_by_user
    default_values["layer_switch_allowed"] = my_filter.layer_switch_allowed_by_user
    default_values["bounds"] = my_filter.bounds_by_user
    #     "calculation_type": my_filter.filter_definition["calculation_type"],
    #     "polarization": my_filter.filter_definition["polarization"],
    #     "polarAngleMin": my_filter.filter_definition["polarAngleMin"],
    #     "polarAngleMax": my_filter.filter_definition["polarAngleMax"],
    #     "polarAngleStep": my_filter.filter_definition["polarAngleStep"],
    #     "azimAngleMin": my_filter.filter_definition["azimAngleMin"],
    #     "azimAngleMax": my_filter.filter_definition["azimAngleMax"],
    #     "azimAngleStep": my_filter.filter_definition["azimAngleStep"],
    #     "wavelengthMin": my_filter.filter_definition["wavelengthMin"],
    #     "wavelengthMax": my_filter.filter_definition["wavelengthMax"],
    #     "wavelengthStep": my_filter.filter_definition["wavelengthStep"],
    #         "structure_materials": my_filter.filter_definition["structure_materials"],
    # "structure_thicknesses": my_filter.filter_definition["structure_thicknesses"],
    # "thickness_opt_allowed": my_filter.filter_definition["thickness_opt_allowed"],
    # "layer_switch_allowed": my_filter.filter_definition["layer_switch_allowed"],
    # "substrate_material": my_filter.filter_definition["substrate_material"],
    # "substrate_thickness": my_filter.filter_definition["substrate_thickness"],
    # "incident_medium": my_filter.filter_definition["incident_medium"],
    # "exit_medium": my_filter.filter_definition["exit_medium"],
    # "targets_type": my_filter.filter_definition["targets_type"],
    # "targets_polarization": my_filter.filter_definition["targets_polarization"],
    # "targets_condition": my_filter.filter_definition["targets_condition"],
    # "targets_value": my_filter.filter_definition["targets_value"],
    # "targets_polar_angle": my_filter.filter_definition["targets_polar_angle"],
    # "polar_angle_steps": my_filter.filter_definition["polar_angle_steps"],
    # "targets_azimuthal_angle": my_filter.filter_definition["targets_azimuthal_angle"],
    # "azimuthal_angle_steps": my_filter.filter_definition["azimuthal_angle_steps"],
    # "targets_wavelengths": my_filter.filter_definition["targets_wavelengths"],
    # "wavelength_steps": my_filter.filter_definition["wavelength_steps"],
    # "targets_tolerance": my_filter.filter_definition["targets_tolerance"],
    #     "bounds": my_filter.filter_definition["bounds"],
    # }

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
        "generalCore": my_filter.filter_definition["core_selection"],
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
    global default_file

    if my_filter is None:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
        ) = upload_file(default_file)
    else:
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            unique_colors,
        ) = extract_filter_design(my_filter)

    return render_template(
        "optimize.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=num_legend_items,
        unique_materials=unique_materials,
        legend_colors=unique_colors,
    )


@app.route("/materials")
def materials():
    # Specify the directory you want to search
    directory = "../materials/"

    # Get a list of all files in the directory
    # file_list = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

    # Get a list of all .csv files in the directory
    material_list = [
        os.path.splitext(f)[0]
        for f in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, f)) and f.endswith(".csv")
    ]

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
        # If provided_filename is provided, use it as the selected file
        selected_file = provided_filename
        filename = provided_filename

    # Create the filter and construct the filter stack representation
    my_filter = FilterStack(selected_file)
    (
        num_boxes,
        colors,
        heights,
        number_unique_materials,
        unique_materials,
        unique_colors,
    ) = extract_filter_design(my_filter)

    if provided_filename is not None:
        return (
            num_boxes,
            colors,
            heights,
            number_unique_materials,
            unique_materials,
            unique_colors,
        )

    else:
        return redirect(url_for("stack"))
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

    ordered_thicknesses = [
        filter.filter_definition["structure_thicknesses"][i] for i in filter.layer_order
    ]

    heights = np.round(np.array(ordered_thicknesses) / sum(ordered_thicknesses) * 300)
    return (
        num_boxes,
        colors.tolist(),
        heights.tolist(),
        len(unique_materials),
        unique_materials.tolist(),
        unique_colors,
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
    data_to_json["substrate_material"] = data[5].get("1").get("values")[0]
    data_to_json["substrate_thickness"] = float(data[6].get("1").get("values")[0])
    data_to_json["incident_medium"] = data[7].get("1").get("values")[0]
    data_to_json["exit_medium"] = data[8].get("1").get("values")[0]
    data_to_json["core_selection"] = (
        "general" if bool(data[9].get("1").get("values")[0]) else "fast"
    )

    layers = np.array(data[10].get("0").get("values"))
    data_to_json["structure_materials"] = layers[0::6].tolist()
    data_to_json["structure_thicknesses"] = np.array(
        [ast.literal_eval(i) if "," in i else float(i) for i in layers[1::6]],
        dtype=object,
    ).tolist()
    data_to_json["thickness_opt_allowed"] = [s.lower() == "true" for s in layers[2::6]]
    data_to_json["bounds"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(layers[3::6].tolist(), layers[4::6].tolist())
    ]
    data_to_json["layer_switch_allowed"] = [s.lower() == "true" for s in layers[5::6]]

    targets = np.array(data[11].get("0").get("values"))
    data_to_json["targets_type"] = targets[0::14].tolist()
    data_to_json["targets_polarization"] = targets[1::14].tolist()
    data_to_json["targets_polar_angle"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(targets[2::14].tolist(), targets[3::14].tolist())
    ]
    data_to_json["polar_angle_steps"] = targets[4::14].astype(float).tolist()
    data_to_json["targets_azimuthal_angle"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(targets[5::14].tolist(), targets[6::14].tolist())
    ]
    data_to_json["azimuthal_angle_steps"] = targets[7::14].astype(float).tolist()
    data_to_json["targets_wavelengths"] = [
        [float(x), float(y)] if y != "" else float(x)
        for x, y in zip(targets[8::14].tolist(), targets[9::14].tolist())
    ]
    data_to_json["wavelength_steps"] = targets[10::14].astype(float).tolist()
    data_to_json["targets_condition"] = targets[11::14].tolist()
    data_to_json["targets_value"] = targets[12::14].astype(float).tolist()
    data_to_json["targets_tolerance"] = targets[13::14].astype(float).tolist()

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
    file_ending = request.args.get("fileEnding")
    if file_ending == ".json":
        # Handle .json file ending
        path_to_file = selected_file
        return send_file(path_to_file, as_attachment=True)
    elif file_ending == ".ofp":
        # Handle .ofp file ending
        # First convert the .json file to .ofp
        input_file = selected_file
        with open(input_file, "r") as input_file:
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
    polarization = data["polarization"]

    my_filter.calculate_ar_data(
        wavelengths,
        polar_angles,
        azimuthal_angles = [float(data["azimuthalAngle"])], 
        target_type=target_type,
        polarization=polarization,
        is_general_core=data["generalCore"],
        web=True,
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
    plot_data = {"x": x, "y": y, "name": str(data["x"])}

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

    merit_time_series = []
    iteration_value = []

    time.sleep(1)

    while not my_filter.stop_flag:
        merit_time_series.append(my_filter.optimum_merit)
        iteration_value.append(my_filter.optimum_iteration)

        # Create a Plotly figure using the data
        # fig = go.Figure(data=go.Scatter(x=iteration_value, y=merit_time_series))
        # Convert the first column to a list and all following columns to a list of lists
        plotting_data = {
            "x": iteration_value,
            "y": merit_time_series,
        }

        # Emit the figure to the client
        socketio.emit("update_merit_graph", plotting_data)

        # Update the stack representation
        (
            num_boxes,
            colors,
            heights,
            number_unique_materials,
            unique_materials,
            unique_colors,
        ) = extract_filter_design(my_filter)

        # Package the values into a dictionary
        filter_representation = {
            "num_boxes": num_boxes,
            "colors": colors,
            "heights": heights,
            "number_unique_materials": number_unique_materials,
            "unique_materials": unique_materials,
            "unique_colors": unique_colors,
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

        time.sleep(1)


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
