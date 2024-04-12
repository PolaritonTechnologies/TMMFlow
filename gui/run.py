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
from flask_socketio import SocketIO, emit
from flask_cors import CORS

# Plotting graphs
import plotly.graph_objects as go
import plotly.io as pio

import threading

# import subprocess
from werkzeug.utils import secure_filename

import os
import numpy as np
import pandas as pd

from utility import (
    allowed_file,
    generate_colors,
)
from FilterStack import FilterStack

from queue import Queue

# Create global queues
message_queue = Queue()
design_update_queue = Queue()
plot_queue = Queue()
plot_done_queue = Queue()

# Create a lock for the queue
message_queue_lock = threading.Lock()
design_update_queue_lock = threading.Lock()
plot_queue_lock = threading.Lock()
plot_done_queue_lock = threading.Lock()


def queue_update():
    print("Queueing design update...")
    with design_update_queue_lock:
        design_update_queue.put(1)


def queue_msg(msg):
    with message_queue_lock:
        message_queue.put(msg)


def queue_plot(plot_dic):
    with plot_queue_lock:
        plot_queue.put(plot_dic)


def queue_plot_done():
    with plot_done_queue_lock:
        plot_done_queue.put(1)


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

num_boxes = None
colors = None
heights = None
unique_materials = None


##############################################
################## Routing ###################
##############################################


@app.route("/", methods=["GET", "POST"])
def home():
    global num_boxes
    global colors
    global heights
    global unique_materials
    global selected_file

    num_boxes, colors, heights, num_legend_items, unique_materials, legend_colors = (
        upload_file(default_file)
    )
    selected_file = default_file
    return render_template(
        "home.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=len(unique_materials),
        unique_materials=unique_materials,
        legend_colors=legend_colors,
        file_label="Currently loaded "
        + default_file
        + ": Click Browse to select another file",
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
    global selected_file
    num_boxes, colors, heights, num_legend_items, unique_materials, legend_colors = (
        upload_file(selected_file)
    )
    return render_template(
        "optimize.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=len(unique_materials),
        unique_materials=unique_materials,
        legend_colors=legend_colors,
    )


##############################################
################## Home ######################
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

    # Create the filter and construct the filter stack representation
    my_filter = FilterStack(selected_file)

    # Now extract the number of layers to construct the number of boxes
    num_boxes = len(my_filter.filter_definition["structure_thicknesses"])
    unique_materials = np.unique(my_filter.filter_definition["structure_materials"])
    unique_colors = generate_colors(len(unique_materials))
    colors = np.empty(num_boxes, dtype=np.dtype("U7"))

    for i in range(len(unique_materials)):
        colors[
            np.array(my_filter.filter_definition["structure_materials"])
            == np.unique(my_filter.filter_definition["structure_materials"])[i]
        ] = unique_colors[i]

    heights = np.round(
        np.array(my_filter.filter_definition["structure_thicknesses"])
        / sum(my_filter.filter_definition["structure_thicknesses"])
        * 400
    )

    if defaultname is not None:
        return (
            num_boxes,
            colors,
            heights,
            len(unique_materials),
            unique_materials,
            unique_colors,
        )

    else:

        return render_template(
            "home.html",
            num_boxes=num_boxes,
            colors=colors,
            heights=heights,
            num_legend_items=len(unique_materials),
            unique_materials=unique_materials,
            legend_colors=unique_colors,
            file_label="Currently loaded: " + filename,
        )


##############################################
############# Calculate and Plot #############
##############################################


@socketio.on("calculate_and_plot")
def calculate_and_plot(data):
    """ """
    global my_filter
    polar_angles = np.arange(float(data["startAngle"]), float(data["endAngle"]) + float(data["stepAngle"]), float(data["stepAngle"]))

    calculated_data_df = pd.DataFrame()

    for theta in polar_angles:
        calculated_data_df[theta]= my_filter.calculate_one_angle(
            float(data["startWavelength"]),
            float(data["endWavelength"]),
            float(data["stepWavelength"]),
            data["mode"],
            data["polarization"],
            theta,
            float(data["azimuthalAngle"]),
            True if data["generalCore"] == "on" else False
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
            colorscale='Viridis',
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
        colorscale='Viridis',
    )

    fig = go.Figure(data=heatmap)

    # Convert the figure to JSON format
    fig_json = pio.to_json(fig)

    # Emit the figure in JSON format
    socketio.emit("update_plot", fig_json)
    


"""
@app.route("/get_plot", methods=["GET"])
def get_plot():
    with plot_queue_lock:
        if plot_queue.empty():
            return jsonify({})
        dict_plot = plot_queue.get()
        print("Updating plot...")
        return jsonify(
            {
                "intensity": dict_plot["intensity"],
                "wavelength": dict_plot["wavelength"],
                "angles": dict_plot["angles"],
                "azimuthal_angle": dict_plot["azimuthal_angle"],
            }
        )
"""


"""
@app.route("/get_plot_done", methods=["GET"])
def get_plot_done():
    with plot_done_queue_lock:
        if plot_done_queue.empty():
            return jsonify({})
        plot_done_queue.get()
        print("Plotting done.")
        return jsonify({"done": True})

"""

##############################################
################ Optimization  ###############
##############################################

@app.route("/get_messages", methods=["GET"])
def get_messages():
    messages = []
    while not message_queue.empty():
        messages.append(message_queue.get())
    return jsonify(messages)


@app.route("/get_design_data", methods=["GET"])
def get_design_data():
    with design_update_queue_lock:
        if design_update_queue.empty():
            return jsonify({})
        design_update_queue.get()
        (
            num_boxes,
            colors,
            heights,
            num_legend_items,
            unique_materials,
            legend_colors,
        ) = upload_file(default_file)
        queue_msg("Updating design...")
        return jsonify(
            {
                "num_boxes": num_boxes,
                "colors": list(colors),
                "heights": list(heights),
                "num_legend_items": len(unique_materials),
                "unique_materials": list(unique_materials),
                "legend_colors": list(legend_colors),
            }
        )


@app.route("/start_optimization", methods=["POST"])
def start_optimization():
    queue_msg("Starting optimization...")
    global stop_optimization
    stop_optimization = False

    # Get the value of optimizationMethod outside of the new thread
    data = request.get_json()
    optimization_method = data.get("optimizationMethod")
    print("Optimization method: ", optimization_method)

    def run_optimization(optimization_method):
        global my_filter
        optim_module = FilterStack(
            "./temp_cpp_order.json",
            my_filter,
            message_queue=message_queue,
            update_queue=design_update_queue,
            log_func=queue_msg,
            log_design_func=queue_update,
        )
        optim_module.perform_optimisation(
            optimization_method,
            save_optimized_to_file=True,
            stop_flag=lambda: stop_optimization,
        )

        # As we are dealing with global variables the new optimized filter has
        # to be specifically assigned to the global variable
        my_filter = optim_module.my_filter

    # Run the optimization in a separate thread
    threading.Thread(target=run_optimization, args=(optimization_method,)).start()

    return "", 204  # Return no content


@app.route("/stop_optimization", methods=["POST"])
def stop_optimization():
    queue_msg("Stopping optimization...")
    global stop_optimization
    stop_optimization = True

    return "", 204  # Return no content


if __name__ == "__main__":
    socketio.run(app)
