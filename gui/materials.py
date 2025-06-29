# Flask and GUI related imports
from flask import (
    Blueprint,
    request,
    jsonify,
    send_file,
    render_template,
    current_app as app,
    session,
    abort,
)
from werkzeug.utils import secure_filename

# basic python packages
import os
import pandas as pd
import numpy as np
import json

# GUI python functions
from .auth import login_required

# from .utility import get_available_materials
from .database import (
    get_available_materials,
    get_material_data,
    deactivate_material_from_db,
    add_material_to_db,
)

# Define the Blueprint
materials_bp = Blueprint("materials_bp", __name__, template_folder="templates")


@materials_bp.route("/materials")
@login_required
def materials():
    """
    Basic view for materials
    """
    # Get a list of all .csv files in the directory
    material_list, material_classes = get_available_materials(session["team"])

    # Group materials by their classes
    grouped_materials = {}
    for material, material_class in zip(material_list, material_classes):
        if material_class not in grouped_materials:
            grouped_materials[material_class] = []
        grouped_materials[material_class].append(material)

    return render_template(
        "materials.html",
        grouped_materials=grouped_materials,
    )


@materials_bp.route("/get_material_data_ajax", methods=["POST"])
def get_material_data_ajax():
    data = request.json
    material = data["material"]

    # Get material data from database and convert back from json to pandas dataframe
    df = get_material_data(material, session["team"])

    # Convert the first column to a list and all following columns to a list of lists
    response_data = {
        "x": df.iloc[:, 0].tolist(),
        "y": df.iloc[:, 1:].values.T.tolist(),
        "name": df.columns[1:].tolist(),
    }

    return jsonify(response_data)


@materials_bp.route("/download_material")
def download_material():
    file_name = request.args.get("fileName")

    # Get material data from database and convert back from json to pandas dataframe
    df = get_material_data(file_name, session["team"])

    # Save the dataframe to a .csv file in the instance/temp folder
    upload_folder = os.path.join(
        os.path.dirname(app.instance_path), app.config["UPLOAD_FOLDER"]
    )
    path_to_file = os.path.join(upload_folder, file_name + ".csv")

    df.to_csv(
        path_to_file,
        sep="\t",
        index=False,
    )

    return send_file(path_to_file, as_attachment=True)


@materials_bp.route("/delete_material", methods=["POST"])
def delete_material():
    material = request.form.get("material")

    deactivate_material_from_db(material, session["team"])

    # Return a success message
    return jsonify({"message": f"Material {material} deleted successfully"})


@materials_bp.route("/upload_material", methods=["POST"])
def upload_material():
    if "file" not in request.files:
        return abort(400, "No file part")
    file = request.files["file"]
    if file.filename == "":
        return abort(400, "No selected file")
    if file:
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config["UPLOAD_FOLDER"], filename))

        # First check the file format (to allow for different types of .csv files)
        with open(os.path.join(app.config["UPLOAD_FOLDER"], filename), "r") as f:
            # Read lines until the first line with numerical data appears
            skiprows = 0
            while True:
                line = f.readline()
                if not line:
                    return abort(400, "File does not contain numerical data")

                if line[0].isdigit():
                    break

                skiprows += 1

            # We want to allow for different character separators (such as \t,
            # ;, ,)
            data_line = f.readline()
            if "\t" in data_line:
                sep = "\t"
            elif ";" in data_line:
                sep = ";"
            elif "," in data_line:
                sep = ","
            else:
                return abort(400, "Data row is not separated by tab, ; or ,")

        # Load file and get data
        df = pd.read_csv(
            os.path.join(app.config["UPLOAD_FOLDER"], filename),
            skiprows=skiprows,
            sep=sep,
            header=None,
        )

        # Delete columns that only contain NaN values
        df = df.dropna(axis=1, how="all")

        # Delete rows that only contain NaN values
        df = df.dropna(axis=0, how="all")

        # Wavelength column
        wl = df.iloc[:, 0]

        # Data columns
        data_columns = df.iloc[:, 1:].values.T

        if np.shape(data_columns)[0] < 2:
            return abort(
                400,
                "Only two columns could be identified (or there is a mismatch between the number of separators of the header and data lines). The data file must contain at least three columns [wavelength, n, k]",
            )

        if np.any(wl > 10000) or np.any(wl < 50):
            return abort(
                400,
                "Wavelengths must be in nm between 50 nm and 10 000 nm. Please check your file format.",
            )

        # Convert the first column to a list and all following columns to a list of lists
        if len(df.columns) == 3:
            name = ["n", "k"]
        elif len(df.columns) == 5:
            name = ["n_ord", "k_ord", "n_exord", "k_exord"]
        elif len(df.columns) == 7:
            name = ["n_x", "k_x", "n_y", "k_y", "n_z", "k_z"]
        else:
            return abort(
                400,
                "The number of columns in the file should be 3 (wavelength, n, k), 5 (wavelenth, n_ord, k_ord, n_exord, k_exord) or 7 (wavelength, n_x, k_x, n_y, k_y, n_z, k_z)",
            )

        data = {
            "x": wl.tolist(),
            "y": data_columns.tolist(),
            "name": name,
        }
        return jsonify({"status": "success", "data": data, "filename": filename}), 200


@materials_bp.route("/accept_material", methods=["POST"])
def accept_material():
    # Read the file in and get the data
    # filename = request.get_json()["filename"]
    plot_data = request.get_json()["plotData"]

    if np.shape(plot_data)[0] == 2:
        df = pd.DataFrame(
            {
                "wavelength": plot_data[0]["x"],
                "n": plot_data[0]["y"],
                "k": plot_data[1]["y"],
            }
        )
    elif np.shape(plot_data)[0] == 4:
        df = pd.DataFrame(
            {
                "wavelength": plot_data[0]["x"],
                "n_ord": plot_data[0]["y"],
                "k_ord": plot_data[1]["y"],
                "n_exord": plot_data[2]["y"],
                "k_exord": plot_data[3]["y"],
            }
        )
    elif np.shape(plot_data)[0] == 6:
        df = pd.DataFrame(
            {
                "wavelength": plot_data[0]["x"],
                "n_x": plot_data[0]["y"],
                "k_x": plot_data[1]["y"],
                "n_y": plot_data[2]["y"],
                "k_y": plot_data[3]["y"],
                "n_z": plot_data[4]["y"],
                "k_z": plot_data[5]["y"],
            }
        )
    else:
        return abort(400, "Invalid number of columns")

    # df = pd.read_csv(
    # os.path.join(app.config["UPLOAD_FOLDER"], filename), skiprows=1, sep="\t"
    # )
    df_as_json = {}
    for col in df.columns:
        df_as_json[col] = df[col].tolist()

    add_material_to_db(
        request.get_json()["materialName"],
        session["team"],
        request.get_json()["materialClass"],
        json.dumps(df_as_json),
        session["user_id"],
    )

    return "File accepted", 200
