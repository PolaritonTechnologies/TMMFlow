from flask import (
    Flask,
    request,
    redirect,
    session,
    render_template,
    flash,
    jsonify,
    url_for,
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
import time

import numpy as np
import pandas as pd
import json

from utility import (
    allowed_file,
    generate_colors,
)
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

    default_values = {
        "num_boxes": num_boxes,
        "colors": colors,
        "heights": heights,
        "num_legend_items": num_legend_items,
        "unique_materials": unique_materials,
        "legend_colors": unique_colors,
        "file_label": default_file,
        "substrate_material_type": my_filter.filter_definition["substrate_material"],
    }


    return render_template(
        "stack.html",
        default_values = default_values,
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

@app.route("/stack_editor")
def stack_editor():
    global my_filter

    if np.all(my_filter.stored_data == None):
        dataPresent = False
    else:
        dataPresent = True

    # default_values = {
    #     "mode": my_filter.filter_definition["calculation_type"],
    #     "startAngle": my_filter.filter_definition["polarAngleMin"],
    #     "endAngle": my_filter.filter_definition["polarAngleMax"],
    #     "stepAngle": my_filter.filter_definition["polarAngleStep"],
    #     "startWavelength": my_filter.filter_definition["wavelengthMin"],
    #     "endWavelength": my_filter.filter_definition["wavelengthMax"],
    #     "stepWavelength": my_filter.filter_definition["wavelengthStep"],
    #     "polarization": my_filter.filter_definition["polarization"],
    #     "azimuthalAngle": my_filter.filter_definition["azimAngleMin"],
    #     "generalCore": my_filter.filter_definition["core_selection"],
    #     "dataPresent": dataPresent,
    # }

    # Send the image URL back to the client
    return render_template(
        "stack_editor.html",
        # default_values=default_values,
    )


##############################################
################## Stack ######################
##############################################


@app.route("/upload_file", methods=["POST"])
def upload_file(defaultname=None):
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
    if defaultname == None:
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
            selected_file = full_path
    else:
        # If defaultname is provided, use it as the selected file
        selected_file = defaultname
        filename = defaultname

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

    if defaultname is not None:
        return (
            num_boxes,
            colors,
            heights,
            number_unique_materials,
            unique_materials,
            unique_colors,
        )

    else:
        return render_template(
            "stack.html",
            num_boxes=num_boxes,
            colors=colors,
            heights=heights,
            num_legend_items=number_unique_materials,
            unique_materials=unique_materials,
            legend_colors=unique_colors,
            file_label="Currently loaded: " + filename,
        )


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

    heights = np.round(
        np.array(filter.filter_definition["structure_thicknesses"])
        / sum(filter.filter_definition["structure_thicknesses"])
        * 300
    )
    return (
        num_boxes,
        colors.tolist(),
        heights.tolist(),
        len(unique_materials),
        unique_materials.tolist(),
        unique_colors,
    )


##############################################
############# Calculate and Plot #############
##############################################


@socketio.on("calculate_and_plot")
def calculate_and_plot(data):
    """ """
    global my_filter
    polar_angles = np.arange(
        float(data["startAngle"]),
        float(data["endAngle"]) + float(data["stepAngle"]),
        float(data["stepAngle"]),
    )

    calculated_data_df = pd.DataFrame()

    for theta in polar_angles:
        calculated_data_df[theta] = my_filter.calculate_one_angle(
            float(data["startWavelength"]),
            float(data["endWavelength"]),
            float(data["stepWavelength"]),
            data["mode"],
            data["polarization"],
            theta,
            float(data["azimuthalAngle"]),
            True if data["generalCore"] == "on" else False,
        )

        # Create a Plotly figure using the calculated data
        angles = calculated_data_df.columns.to_numpy()
        wavelengths = calculated_data_df.index.to_numpy()
        color_values = calculated_data_df.to_numpy()

        # The layout is stored in simulate.html
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

    my_filter.stored_data = calculated_data_df


@socketio.on("plot")
def plot():
    """ """
    global my_filter

    calculated_data_df = my_filter.stored_data

    # Create a Plotly figure using the calculated data
    angles = calculated_data_df.columns.to_numpy()
    wavelengths = calculated_data_df.index.to_numpy()
    color_values = calculated_data_df.to_numpy()

    # The layout is stored in simulate.html
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
        fig = go.Figure(data=go.Scatter(x=iteration_value, y=merit_time_series))

        # Convert the figure to JSON format
        fig_json = pio.to_json(fig)

        # Emit the figure to the client
        socketio.emit("update_merit_graph", fig_json)

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
            'num_boxes': num_boxes,
            'colors': colors,
            'heights': heights,
            'number_unique_materials': number_unique_materials,
            'unique_materials': unique_materials,
            'unique_colors': unique_colors
        }

        # Convert the dictionary to a JSON string
        filter_json = json.dumps(filter_representation)

        # Emit the JSON string
        socketio.emit("update_filter_representation", filter_json)

        time.sleep(1)


@socketio.on("stop_optimization")
def stop_optimization():
    """ """
    global my_filter

    my_filter.stop_flag = True


if __name__ == "__main__":
    socketio.run(app)
