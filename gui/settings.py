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
import pandas as pd

# GUI python functions
from .auth import login_required
from .database import change_password_database

# Define the Blueprint
settings_bp = Blueprint("settings_bp", __name__, template_folder="templates")


@settings_bp.route("/settings")
@login_required
def settings():
    """
    !# should be updated so that the materials are also read in from a database at some point
    Basic view for materials
    """
    # Get a list of all .csv files in the directory
    default_values = {
        "team": session["team"],
        "email": session["email"],
        "username": session["user_id"],
    }

    return render_template("settings.html", default_values=default_values)


@settings_bp.route("/change_password", methods=["POST"])
def change_password():
    password = request.json["passwordChanged"]
    password_repeated = request.json["passwordRepeated"]

    if password != password_repeated:
        return jsonify({"status": "error", "message": "Passwords do not match"})
    else:
        change_password_database(session["user_id"], password)

        response_data = {
            "status": "success",
            "message": "Password changed successfully",
        }

        return jsonify(response_data)


@settings_bp.route("/download_material")
def download_material():
    file_name = request.args.get("fileName")
    path_to_file = app.config["MATERIAL_FOLDER"] + file_name + ".csv"

    return send_file(path_to_file, as_attachment=True)


"""
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

        # Emit the data
        #! Check if this is still working without socket.io
        # socketio.emit("material_data", data)

        return filename, 200
"""
