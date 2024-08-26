# Flask and GUI related imports
from flask_login import login_required
from flask import (
    Blueprint,
    json,
    request,
    redirect,
    render_template,
    jsonify,
    url_for,
    send_file,
    session,
    current_app as app,
)
from werkzeug.utils import secure_filename

# Calculation and basic python packages
import numpy as np
import os
from datetime import datetime
import ast

# GUI python functions
from .database import (
    Job,
    get_latest_job_id,
    load_latest_filter,
    select_latest_optimization,
    select_initial_optimization,
    get_all_filter_versions,
    select_optimization_by_datetime,
    get_all_user_projects,
    select_job_by_datetime_and_name,
)
from . import db
from .utility import (
    get_available_materials,
    get_available_templates,
    allowed_file,
    generate_colors,
    is_number,
)

# Core python functions
from open_filter_converter import import_from_open_filter, export_to_open_filter
from FilterStack import FilterStack

# Define the Blueprint
stack_bp = Blueprint("stack_bp", __name__, template_folder="templates")


@stack_bp.route("/stack", methods=["GET", "POST"])
@login_required
def stack():
    """
    View for the basic overview over the filter stack design.
    """
    # If no filter has been initialized so far, load the default file
    if session.get("job_id") is None:
        # The upload_file function calls stack at the end of the function again
        # and sets the filter_initialized flag to True
        upload_file(app.config["DEFAULT_FILE"])

    # Read in the latest design of this session
    job_id = session.get("job_id")
    filter_definition = load_latest_filter(job_id).return_current_design_as_json()
    (
        num_boxes,
        colors,
        heights,
        num_legend_items,
        unique_materials,
        unique_colors,
        incoherent,
    ) = extract_filter_design(
        filter_definition["structure_thicknesses"],
        filter_definition["structure_materials"],
        filter_definition["incoherent"],
    )

    # Read in the latest optimization of this session
    optimization = select_latest_optimization(session.get("job_id"))
    current_structure = json.loads(optimization.current_json)

    # Specify the directory you want to search
    material_list = get_available_materials(app.config["MATERIAL_FOLDER"])
    template_list = get_available_templates(app.config["TEMPLATE_FOLDER"])

    # Get all versions of current filter design
    filter_versions = get_all_filter_versions(job_id)
    filter_versions_time_stamps = [
        filter_version.time_stamp.strftime("%Y-%m-%d %H:%M:%S")
        for filter_version in filter_versions
    ]

    # Read in all available projects in the past for each user
    unique_job_ids, unique_filter_names, unique_timestamps = get_all_user_projects(
        session.get("user_id")
    )

    default_values = {
        "num_boxes": num_boxes,
        "colors": colors,
        "heights": heights,
        "num_legend_items": num_legend_items,
        "unique_materials": unique_materials,
        "incoherent_boxes": incoherent,
        "legend_colors": unique_colors,
        "file_label": session.get("filter_name"),
        "available_materials": material_list,
        "available_templates": template_list,
        "versions": filter_versions_time_stamps,
        "previous_projects": [
            unique_timestamps[i].strftime("%Y-%m-%d %H:%M:%S")
            + ": "
            + unique_filter_names[i]
            for i in range(len(unique_job_ids))
        ],
    }

    # Add the entire filter definition to the default values to render it
    default_values.update(current_structure)
    temp_material_list = current_structure["structure_materials"]
    default_values["structure_materials"] = [
        (
            item.split("_")[0] + "_" + "_".join(item.split("_")[1:][::-1])
            if "_" in item
            else item
        )
        for item in temp_material_list
    ]
    default_values["structure_thicknesses"] = [
        np.flip(a).tolist() for a in current_structure["structure_thicknesses"]
    ]
    default_values["thickness_opt_allowed"] = current_structure["thickness_opt_allowed"]
    default_values["layer_switch_allowed"] = current_structure["layer_switch_allowed"]
    default_values["bounds"] = current_structure["bounds"]
    default_values["incoherent"] = current_structure["incoherent"]

    return render_template(
        "stack.html",
        default_values=default_values,
    )


def check_uploaded_filter(filter_definition_json):
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
        return dict(
            {
                "error": "Some materials are not available",
                "mismatched_materials": unavailable_materials,
                "available_materials": available_materials,
            }
        )
    else:
        # Now return no error and the updated filter definition
        return dict(
            {"error": None, "updated_filter_definition": filter_definition_json}
        )


@stack_bp.route("/upload_file", methods=["POST"])
def upload_file(filename=None, filter_name=None):
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
    my_filter_path = None

    # User upload of a file
    # Check if the request contains a file
    if filename == None:
        if "file" in request.files:
            file = request.files["file"]
            checked_filename = secure_filename(file.filename)
            full_path = os.path.join(app.config["UPLOAD_FOLDER"], checked_filename)
            file.save(full_path)

        # Check if the checked_filename is allowed
        if allowed_file(checked_filename):
            full_path = os.path.join(app.config["UPLOAD_FOLDER"], checked_filename)
            if checked_filename.split(".")[-1] == "json":
                # Standard .json format
                my_filter_path = full_path
            elif checked_filename.split(".")[-1] == "ofp":
                # Open filter format (has to be translated to json first)
                json_path = import_from_open_filter(full_path)
                my_filter_path = json_path

    else:
        # If the request type is a get, use the default file instead that was
        # passed as argument
        my_filter_path = filename

    # Now check if the selected file is valid
    with open(my_filter_path, "r") as json_file_path:
        filter_definition_json = json.load(json_file_path)

    response = check_uploaded_filter(filter_definition_json)

    if response["error"] is not None:
        return jsonify(response)
    else:
        filter_definition_json = response["updated_filter_definition"]

    # Create a new job corresponding to a new filter in the job database
    if filter_name == None:
        session["filter_name"] = my_filter_path.split("/")[-1].split(".")[0]
    else:
        session["filter_name"] = filter_name

    # Also create a new optimization corresponding to the new filter in the
    # optimization database with optimization method "None" to track the initial
    # state
    new_job = Job(
        job_id=get_latest_job_id() + 1,
        username=session["user_id"],
        filter_name=session["filter_name"],
        time_stamp=datetime.now(),
        optimization_method="None",
        current_json=json.dumps(filter_definition_json),
    )
    db.session.add(new_job)
    db.session.commit()

    # store the job_id inside the session. It is the unique identifier for the
    # user session and used to retrieve data from the database throughout the
    # software. Everytime a new filter is loaded, a new job_id is created.
    session["job_id"] = new_job.job_id

    # If upload was called from a python function return nothing, else redirect to stack
    if filename is None:
        return jsonify({"redirect": url_for("stack_bp.stack")})
    else:
        return


def extract_filter_design(thicknesses, materials, incoherent_list):
    """
    Extracts the filter design information.

    Args:
        filter_json: the json or path to the json.

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
    num_boxes = len(thicknesses)
    unique_materials = np.unique(materials)
    unique_colors = generate_colors(len(unique_materials))
    colors = np.empty(num_boxes, dtype=np.dtype("U7"))

    for i in range(len(unique_materials)):
        colors[np.array(materials) == np.unique(materials)[i]] = unique_colors[i]

    # Set the incoherent layer thicknesses to a specific value so that it does
    # not go out of control thick
    thicknesses = np.array(thicknesses)
    thicknesses[np.array(incoherent_list)] = 100

    heights = np.round(np.array(thicknesses) / sum(thicknesses) * 300)
    return (
        num_boxes,
        colors.tolist(),
        heights.tolist(),
        len(unique_materials),
        unique_materials.tolist(),
        unique_colors,
        incoherent_list,
    )


@stack_bp.route("/save_json", methods=["POST"])
def save_json():
    """
    Function to read in the current filter design from the GUI and save it to
    the database
    """

    data = request.get_json()
    # button = request.args.get("button")

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
            [
                ast.literal_eval(i) if "," in i else float(i)
                for i in layers[1::number_of_columns]
            ],
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
            for x, y in zip(
                layers[3::number_of_columns].tolist(),
                layers[4::number_of_columns].tolist(),
            )
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
    if np.size(targets) == 1:
        # In case no targets were defined
        data_to_json["targets_type"] = []
        data_to_json["targets_polarization"] = []
        data_to_json["targets_polar_angle"] = []
        data_to_json["polar_angle_steps"] = []
        data_to_json["targets_azimuthal_angle"] = []
        data_to_json["azimuthal_angle_steps"] = []
        data_to_json["targets_wavelengths"] = []
        data_to_json["wavelength_steps"] = []
        data_to_json["targets_condition"] = []
        data_to_json["targets_value"] = []
        data_to_json["targets_tolerance"] = []
        data_to_json["targets_arithmetic"] = []
    else:
        data_to_json["targets_type"] = targets[0::no_columns_targets].tolist()
        data_to_json["targets_polarization"] = targets[1::no_columns_targets].tolist()
        data_to_json["targets_polar_angle"] = [
            [float(x), float(y)] if y != "" else float(x)
            for x, y in zip(
                targets[2::no_columns_targets].tolist(),
                targets[3::no_columns_targets].tolist(),
            )
        ]
        data_to_json["polar_angle_steps"] = (
            targets[4::no_columns_targets].astype(float).tolist()
        )
        data_to_json["targets_azimuthal_angle"] = [
            [float(x), float(y)] if y != "" else float(x)
            for x, y in zip(
                targets[5::no_columns_targets].tolist(),
                targets[6::no_columns_targets].tolist(),
            )
        ]
        data_to_json["azimuthal_angle_steps"] = (
            targets[7::no_columns_targets].astype(float).tolist()
        )
        data_to_json["targets_wavelengths"] = [
            [float(x), float(y)] if y != "" else float(x)
            for x, y in zip(
                targets[8::no_columns_targets].tolist(),
                targets[9::no_columns_targets].tolist(),
            )
        ]
        data_to_json["wavelength_steps"] = (
            targets[10::no_columns_targets].astype(float).tolist()
        )
        data_to_json["targets_condition"] = targets[11::no_columns_targets].tolist()
        data_to_json["targets_value"] = (
            targets[12::no_columns_targets].astype(float).tolist()
        )
        data_to_json["targets_tolerance"] = (
            targets[13::no_columns_targets].astype(float).tolist()
        )
        data_to_json["targets_arithmetic"] = (
            targets[14::no_columns_targets].astype(str).tolist()
        )

    # Now save the data to the database (in the optimization database)
    optimization = Job(
        job_id=session.get("job_id"),
        username=session.get("user_id"),
        filter_name=session.get("filter_name"),
        optimization_method="None",
        current_json=json.dumps(data_to_json),
        time_stamp=datetime.now(),
    )
    db.session.add(optimization)
    db.session.commit()

    # Reload the page
    stack()
    return "", 302


@stack_bp.route("/download")
def download_file():
    """
    Function to download the current
    """

    file_type = request.args.get("fileType")

    if file_type == "Initial":
        optimization = select_initial_optimization(session.get("job_id"))
        json_structure = json.loads(optimization.current_json)
    elif file_type == "Optimized":
        optimization = select_latest_optimization(session.get("job_id"))
        json_structure = json.loads(optimization.current_json)

    # Dump it to a json file
    file_ending = request.args.get("fileEnding")
    path_to_file = os.path.join(
        os.getcwd(),
        app.config["UPLOAD_FOLDER"],
        file_type + "_" + str(session.get("job_id")),
    )

    if file_ending == ".json":
        # Handle .json file ending
        # Dump json into a file
        with open(path_to_file + file_ending, "w") as f:
            json.dump(json_structure, f)
        return send_file(path_to_file + file_ending, as_attachment=True)
    elif file_ending == ".ofp":
        # Handle .ofp file ending
        # First convert the cpp .json file to .ofp
        path_to_file = export_to_open_filter(
            FilterStack(json_structure).return_current_design_as_json(),
            path_to_file.split("/")[-1],
        )
        return send_file(path_to_file, as_attachment=True)


@stack_bp.route("/reset_filter", methods=["POST"])
def reset_filter():
    # Copy the first database entry for the given key to a new optimization entry (the newest)
    optimization = select_initial_optimization(session.get("job_id"))
    initial_json = json.loads(optimization.current_json)

    # Now create a new optimization entry with the initial json
    new_optimization = Job(
        job_id=session.get("job_id"),
        username=session.get("user_id"),
        filter_name=session.get("filter_name"),
        optimization_method="None",
        current_json=json.dumps(initial_json),
        time_stamp=datetime.now(),
    )

    # Commit the new optimization to the database
    db.session.add(new_optimization)
    db.session.commit()

    # Reload page
    return redirect("/")


@stack_bp.route("/select_version", methods=["POST"])
def select_version():
    selected_version = datetime.strptime(
        request.form.get("version"), "%Y-%m-%d %H:%M:%S"
    )

    # Retrieve the optimization entry for the selected version
    optimization = select_optimization_by_datetime(
        selected_version, session.get("job_id")
    )
    initial_json = json.loads(optimization.current_json)

    # Now create a new optimization entry with the initial json
    new_optimization = Job(
        job_id=session.get("job_id"),
        username=session.get("user_id"),
        filter_name=session.get("filter_name"),
        optimization_method="None",
        current_json=json.dumps(initial_json),
        time_stamp=datetime.now(),
    )

    # Commit the new optimization to the database
    db.session.add(new_optimization)
    db.session.commit()

    # Reload page
    return redirect("/")


@stack_bp.route("/start_new_design", methods=["POST"])
def start_new_design():
    # Load in template
    template = request.form.get("template") + ".json"
    filter_name = request.form.get("filter_name")

    # Trigger upload file function with the selected template
    upload_file(os.path.join(app.config["TEMPLATE_FOLDER"], template), filter_name)

    return jsonify({"redirect": url_for("stack_bp.stack")})


@stack_bp.route("/load_design", methods=["POST"])
def load_design():
    # Load in template
    filter_version = request.form.get("filter_version")

    # Convert the select elements text to a datetime object and a string with
    # the filter name
    filter_date = filter_version.split(" ")[0] + " " + filter_version.split(" ")[1]
    filter_date = datetime.strptime(filter_date[:-1], "%Y-%m-%d %H:%M:%S")
    filter_name = filter_version.split(" ")[-1]

    # Retrieve the optimization entry for the selected version
    job_state = select_job_by_datetime_and_name(filter_date, filter_name)

    # Now set the job_id to what was selected and reload the page
    session["job_id"] = job_state.job_id
    session["filter_name"] = filter_name

    return jsonify({"redirect": url_for("stack_bp.stack")})
