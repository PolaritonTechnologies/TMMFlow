# Flask and GUI related imports
from flask import (
    Blueprint,
    request,
    jsonify,
    send_file,
    render_template,
    current_app as app,
    session,
)
from werkzeug.utils import secure_filename

# basic python packages
import os
import shutil
import pandas as pd

# GUI python functions
from .auth import login_required

# from .utility import get_available_materials
from .database import (
    get_available_materials,
    get_material_data,
    delete_material_from_db,
)

# Define the Blueprint
materials_bp = Blueprint("materials_bp", __name__, template_folder="templates")


@materials_bp.route("/materials")
@login_required
def materials():
    """
    !# should be updated so that the materials are also read in from a database at some point
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

    delete_material_from_db(material, session["team"])

    # Return a success message
    return jsonify({"message": f"Material {material} deleted successfully"})


@materials_bp.route("/upload_material", methods=["POST"])
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
        return jsonify({"status": "success", "data": data, "filename": filename}), 200


@materials_bp.route("/accept", methods=["POST"])
def accept_file():
    filename = request.form.get("filename")
    shutil.move(
        os.path.join(app.config["UPLOAD_FOLDER"], filename),
        os.path.join(app.config["MATERIAL_FOLDER"], filename),
    )
    return "File accepted", 200


@materials_bp.route("/reject", methods=["POST"])
def reject_file():
    filename = request.form.get("filename")
    os.remove(os.path.join(app.config["UPLOAD_FOLDER"], filename))
    return "File rejected", 200
