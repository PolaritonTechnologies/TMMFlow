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

# Create a lock for the queue
design_update_queue_lock = threading.Lock()


def queue_update(msg):
    with design_update_queue_lock:

        if not design_update_queue.empty():
            design_update_queue.put(1)


def queue_msg(msg):
    message_queue.put(msg)


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


@app.route("/", methods=["GET", "POST"])
def home():
    global num_boxes
    global colors
    global heights
    global unique_materials

    num_boxes, colors, heights, num_legend_items, unique_materials, legend_colors = (
        upload_file(default_file)
    )
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
    return render_template("simulate.html")


@app.route("/optimize", methods=["GET", "POST"])
def optimize():
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
    # Perform calculations here...
    plotting = CalculationModule(my_filter, lib)

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

    print("plotting results...")
    plotting.calculate_ar_data(
        wavelength,
        polar_angles,
        azimuthal_angles,
        "r",
        "s",
        save_figure=True,
        save_data=True,
    )

    # Send the image URL back to the client
    return render_template("simulate.html", plot_url="plot.png")


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
        else:
            design_update_queue.get()
            (
                num_boxes,
                colors,
                heights,
                num_legend_items,
                unique_materials,
                legend_colors,
            ) = upload_file("current_structure.json")
        return jsonify(
            {
                "num_boxes": num_boxes,
                "colors": colors,
                "heights": heights,
                "num_legend_items": len(unique_materials),
                "unique_materials": unique_materials,
                "legend_colors": colors,
            }
        )


@app.route("/start_optimization", methods=["POST"])
def start_optimization():
    queue_msg("Starting optimization...")
    global stop_optimization
    stop_optimization = False

    # Get the value of optimizationMethod outside of the new thread
    optimization_method = request.form.get("optimizationMethod")
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

    return render_template(
        "optimize.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=len(unique_materials),
        unique_materials=unique_materials,
        legend_colors=colors,
        optimization_running=True,
    )


@app.route("/stop_optimization", methods=["POST"])
def stop_optimization():
    queue_msg("Stopping optimization...")
    global stop_optimization
    stop_optimization = True
    return render_template(
        "optimize.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=len(unique_materials),
        unique_materials=unique_materials,
        legend_colors=colors,
        optimization_running=False,
    )


if __name__ == "__main__":
    app.run(debug=True)
