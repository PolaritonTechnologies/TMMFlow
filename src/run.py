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
import threading

# import subprocess
from werkzeug.utils import secure_filename

import json
import os
import numpy as np

from utility import (
    translate_order_for_cpp,
    create_filter,
    allowed_file,
    generate_colors,
)
from optim_module import OptimModule
from calculation_module import CalculationModule

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


app = Flask(__name__)
app.secret_key = "your secret key"  # replace with your secret key

my_filter = None
lib = None
optimisation_order_file = None
default_file = "current_structure.json"
num_boxes = None
colors = None
heights = None
unique_materials = None
selected_file = None


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
        legend_colors=colors,
        file_label="Currently loaded "
        + default_file
        + ": Click Browse to select another file",
    )


@app.route("/simulate", methods=["GET", "POST"])
def simulate():

    with open(optimisation_order_file) as f:
        plot_order = json.load(f)

    wavelength = np.arange(
        plot_order["wavelengthMin"],
        plot_order["wavelengthMax"] + 1,
        plot_order["wavelengthStep"],
    )
    polar_angles = np.arange(
        plot_order["polarAngleMin"],
        plot_order["polarAngleMax"] + 1,
        plot_order["polarAngleStep"],
    )
    azimuthal_angles = np.arange(
        plot_order["azimAngleMin"],
        plot_order["azimAngleMax"] + 1,
        plot_order["azimAngleStep"],
    )
    # Send the image URL back to the client
    return render_template(
        "simulate.html",
        intensity=[
            [0 for _ in range(len(polar_angles))] for _ in range(len(wavelength))
        ],
        wavelength=list(wavelength),
        angles=list(polar_angles),
        azimuthal_angle=0,
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
        legend_colors=colors,
    )


@app.route("/upload_file", methods=["POST"])
def upload_file(defaultname=None):
    global selected_file
    if defaultname == None:
        if "file" not in request.files:
            flash("No file part")
            return redirect(request.url)
        file = request.files["file"]
        if file.filename == "":
            flash("No selected file")
            return redirect(request.url)
    else:
        file = open(defaultname, "r")
        file.filename = defaultname

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        selected_file = filename
        data = json.load(file)

        # Save the data back to a JSON file so that it can be accessed by the C++ code
        with open(os.path.join("./", "output.json"), "w") as output_file:
            json.dump(data, output_file)

        global my_filter
        global lib
        global optimisation_order_file

        optimisation_order_file_python = "output.json"
        optimisation_order_file = translate_order_for_cpp(
            optimisation_order_file_python
        )
        my_filter, lib = create_filter(optimisation_order_file)

        # Now extract the number of layers to construct the number of boxes
        with open(os.path.join("./", "temp_cpp_order.json"), "r") as output_file:
            cpp_rendered_data = json.load(output_file)

        num_boxes = len(cpp_rendered_data["structure_thicknesses"])
        unique_materials = np.unique(cpp_rendered_data["structure_materials"])
        unique_colors = generate_colors(len(unique_materials))
        colors = np.empty(num_boxes, dtype=np.dtype("U7"))

        for i in range(len(unique_materials)):
            colors[
                np.array(cpp_rendered_data["structure_materials"])
                == np.unique(cpp_rendered_data["structure_materials"])[i]
            ] = unique_colors[i]

        heights = np.round(
            np.array(cpp_rendered_data["structure_thicknesses"])
            / sum(cpp_rendered_data["structure_thicknesses"])
            * 400
        )

        file.close()

        if defaultname is not None:
            return (
                num_boxes,
                colors,
                heights,
                len(unique_materials),
                unique_materials,
                colors,
            )

        else:

            return render_template(
                "home.html",
                num_boxes=num_boxes,
                colors=colors,
                heights=heights,
                num_legend_items=len(unique_materials),
                unique_materials=unique_materials,
                legend_colors=colors,
                file_label="Currently loaded: " + filename,
            )


@app.route("/calculate_and_plot", methods=["POST"])
def calculate_and_plot():

    plotting = CalculationModule(
        my_filter,
        lib,
        plot_queue,
        plot_done_queue,
        log_plot=queue_plot,
        log_plot_done=queue_plot_done,
    )

    with open(optimisation_order_file) as f:
        plot_order = json.load(f)

    wavelength = np.arange(
        plot_order["wavelengthMin"],
        plot_order["wavelengthMax"] + 1,
        plot_order["wavelengthStep"],
    )
    polar_angles = np.arange(
        plot_order["polarAngleMin"],
        plot_order["polarAngleMax"] + 1,
        plot_order["polarAngleStep"],
    )
    azimuthal_angles = np.arange(
        plot_order["azimAngleMin"],
        plot_order["azimAngleMax"] + 1,
        plot_order["azimAngleStep"],
    )

    plotting.calculate_ar_data(
        wavelength,
        polar_angles,
        azimuthal_angles,
        plot_order["calculation_type"],
        plot_order["polarization"],
        save_figure=True,
        save_data=True,
    )

    return "", 204  # Return no content


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


@app.route("/get_plot_done", methods=["GET"])
def get_plot_done():
    with plot_done_queue_lock:
        if plot_done_queue.empty():
            return jsonify({})
        plot_done_queue.get()
        print("Plotting done.")
        return jsonify({"done": True})


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
        ) = upload_file("current_structure.json")
        queue_msg("Updating design...")
        return jsonify(
            {
                "num_boxes": num_boxes,
                "colors": list(colors),
                "heights": list(heights),
                "num_legend_items": len(unique_materials),
                "unique_materials": list(unique_materials),
                "legend_colors": list(colors),
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
        optim_module = OptimModule(
            "temp_cpp_order.json",
            my_filter,
            lib,
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
    app.run(debug=True)
