# Flask and GUI related imports
from flask import (
    Blueprint,
    json,
    request,
    jsonify,
    send_file,
    session,
    render_template,
    current_app as app,
)

# Calculation and basic python packages
import numpy as np
import pandas as pd
import os

# GUI python functions
from . import db
from .database import load_latest_filter, select_latest_optimization
from .auth import login_required

# Define the Blueprint
simulate_bp = Blueprint("simulate_bp", __name__, template_folder="templates")


@simulate_bp.route("/simulate")
@login_required
def simulate():
    """
    Simulate view
    """
    optimization = select_latest_optimization(session.get("job_id"))
    current_json = json.loads(optimization.current_json)

    if np.all(optimization.current_data is None):
        dataPresent = False
    else:
        dataPresent = True

    default_values = {
        "mode": current_json["calculation_type"],
        "startAngle": current_json["polarAngleMin"],
        "endAngle": current_json["polarAngleMax"],
        "stepAngle": current_json["polarAngleStep"],
        "startWavelength": current_json["wavelengthMin"],
        "endWavelength": current_json["wavelengthMax"],
        "stepWavelength": current_json["wavelengthStep"],
        "polarization": current_json["polarization"],
        "azimuthalAngle": current_json["azimAngleMin"],
        "dataPresent": dataPresent,
    }

    # Send the image URL back to the client
    return render_template(
        "simulate.html",
        default_values=default_values,
    )


@simulate_bp.route("/calculate_and_plot_ajax", methods=["POST"])
def calculate_and_plot_ajax():
    data = request.json

    job_id = session.get("job_id")
    my_filter = load_latest_filter(job_id)
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

    angles = calculated_data_df.columns.to_numpy()
    wavelengths = calculated_data_df.index.to_numpy()
    color_values = calculated_data_df.to_numpy()

    plotting_data = {
        "x": angles.tolist(),
        "y": wavelengths.tolist(),
        "z": color_values.tolist(),
    }

    # Now update the current_data in the optimization database by adding it to
    # the latest optimization entry
    select_latest_optimization(job_id).current_data = json.dumps(plotting_data)
    db.session.commit()

    current_json = my_filter.filter_definition

    additional_data = extract_spectral_parameters(current_json, plotting_data)
    plotting_data.update(additional_data)

    return jsonify(plotting_data)


@simulate_bp.route("/plot_data", methods=["POST"])
def plot_data():
    # Check if there is data available for the job_id in the database
    optimization = select_latest_optimization(session.get("job_id"))
    if optimization.current_data is not None:
        plotting_data = json.loads(optimization.current_data)
        current_json = json.loads(optimization.current_json)

        additional_data = extract_spectral_parameters(current_json, plotting_data)
        plotting_data.update(additional_data)

        # Return the plotting data as a JSON response
        return jsonify(plotting_data)


@simulate_bp.route("/plot_xy_ajax", methods=["POST"])
def plot_xy_ajax():
    # Retrieve the data from the request
    request_data = request.json

    # Load the already calculated data from the sql database using the unique
    # identifier job_id and selecting the latest entry in the optimization table
    optimization = select_latest_optimization(session.get("job_id"))
    plotting_data = json.loads(optimization.current_data)

    # Extract the y data given the point that was clicked on defined by
    # data["x"]
    x = plotting_data["y"]
    y = (
        np.array(plotting_data["z"])
        .T[np.where(np.isclose(plotting_data["x"], request_data["x"]))[0][0]]
        .tolist()
    )

    # x = my_filter.stored_data[0].index.to_list()

    # Get the y data corresponding to the x data
    # y = my_filter.stored_data[0].loc[:, data["x"]].to_list()

    # Create a dictionary with the x and y data
    plot_data = {
        "x": x,
        "y": y,
        "name": str(request_data["x"])
        + "°, "
        + str(request_data["mode"])
        + ", "
        + str(request_data["polarization"]),
    }

    return jsonify(plot_data)


@simulate_bp.route("/download_data")
def download_data():
    optimization = select_latest_optimization(session.get("job_id"))
    plotting_data = json.loads(optimization.current_data)

    # Now dump this into a file with the job_id as the name
    file_path = os.path.join(
        os.getcwd(),
        app.config["UPLOAD_FOLDER"],
        "simulated_data_" + str(session.get("job_id")) + ".csv",
    )
    pd.DataFrame(
        plotting_data["z"], columns=plotting_data["x"], index=plotting_data["y"]
    ).to_csv(file_path, sep="\t")

    return send_file(file_path, as_attachment=True)


def extract_spectral_parameters(current_json, plotting_data):
    target_passband_wavelength = np.array(
        current_json["targets_wavelengths"], dtype=object
    )[np.array(current_json["targets_value"]) > 0.1]

    # Integrate the spectrum for the target passband
    # Only choose first target (target_no = 1)
    target_no = 0
    if np.size(target_passband_wavelength[target_no]) == 1:
        integrated_passband = 0
    else:
        # Check that the target passband is within the range of the plotted data
        if target_passband_wavelength[target_no][1] > np.max(plotting_data["y"]):
            integrated_passband = 0
        elif target_passband_wavelength[target_no][0] < np.min(plotting_data["y"]):
            integrated_passband = 0
        else:
            target_index_low = np.where(
                np.isclose(plotting_data["y"], target_passband_wavelength[target_no][0])
            )[0][0]

            target_index_high = np.where(
                np.isclose(plotting_data["y"], target_passband_wavelength[target_no][1])
            )[0][0]
            integrated_passband = np.trapz(
                [
                    np.trapz(
                        np.array(plotting_data["z"]).T[i][
                            target_index_low:target_index_high
                        ],
                        np.array(
                            plotting_data["y"][target_index_low:target_index_high]
                        ),
                    )
                    * 100
                    for i in range(np.size(plotting_data["x"]))
                ],
                plotting_data["x"],
            )

    # Do calculations on the plotting data
    integrated_spectrum = np.trapz(
        [
            np.trapz(np.array(plotting_data["z"]).T[i], np.array(plotting_data["y"]))
            * 100
            for i in range(np.size(plotting_data["x"]))
        ],
        plotting_data["x"],
    )
    peak_spectrum = np.max(plotting_data["z"]) * 100

    additional_data = {
        "integrated_spectrum": integrated_spectrum,
        "integrated_passband": integrated_passband,
        "peak_spectrum": peak_spectrum,
    }

    return additional_data
