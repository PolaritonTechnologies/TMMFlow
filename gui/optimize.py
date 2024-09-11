# Flask and GUI related imports
from flask import (
    Blueprint,
    json,
    request,
    jsonify,
    session,
    render_template,
)

from datetime import datetime

from . import db, q, redis_conn
from .auth import login_required
from .database import load_latest_filter, select_latest_optimization, Job
from .stack import extract_filter_design

from FilterStack import FilterStack, optimization_function

import rq
import numpy as np

# Define the Blueprint
optimize_bp = Blueprint("optimize_bp", __name__, template_folder="templates")


@optimize_bp.route("/optimize", methods=["GET", "POST"])
@login_required
def optimize():
    """
    Basic view for optimization
    """
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

    # Check if there is already an optimization running
    rq_job_id = redis_conn.get(f"custom_job_id:{job_id}")

    try:
        job = rq.job.Job.fetch(rq_job_id.decode("utf-8"), connection=redis_conn)

        if job.is_started:
            print("A job is currently running.")
            job_status = "running"
        elif job.is_queued:
            print("A job is currently queued.")
            job_status = "queued"
        else:
            print("The job is not running.")
            job_status = "not running"
    except (rq.exceptions.NoSuchJobError, AttributeError):
        print("Currently there is no optimization running for this job!")
        job_status = "not running"

    return render_template(
        "optimize.html",
        num_boxes=num_boxes,
        colors=colors,
        heights=heights,
        num_legend_items=num_legend_items,
        unique_materials=unique_materials,
        legend_colors=unique_colors,
        incoherent=incoherent,
        job_status=job_status,
    )


@optimize_bp.route("/start_optimization", methods=["POST"])
def start_optimization():
    job_id = session.get("job_id")

    # Get the latest optimization entry from the database and create the actual
    # filter object
    latest_design = json.loads(select_latest_optimization(job_id).current_json)

    # Also add a new entry to the optimization table with the optimization
    # method that will be updated throughout the optimization process
    filter = FilterStack(latest_design)
    initial_merit = filter.calculate_initial_merit()
    del filter

    # Set the initial merit and current merit so that the updates are not taking
    # previous values
    redis_conn.set(f"current_:{job_id}", "false")
    intermediate_result = {
        "step": 0,
        "merit": initial_merit,
    }
    redis_conn.set(
        f"current:{job_id}",
        json.dumps(intermediate_result),
    )
    intermediate_result["current_structure"] = latest_design
    redis_conn.set(
        f"current_best:{job_id}",
        json.dumps(intermediate_result),
    )

    optimization = Job(
        job_id=job_id,
        time_stamp=datetime.now(),
        username=session.get("user_id"),
        filter_name=session.get("filter_name"),
        optimization_method=json.dumps(request.json["optimizationMethod"]),
        initial_merit=initial_merit,
        current_json=json.dumps(latest_design),
        description="Optimization with " + str(request.json["optimizationMethod"]),
    )
    db.session.add(optimization)
    db.session.commit()

    # Set the cancellation flag in the redis database to false
    redis_conn.set(f"job_cancel_flag:{job_id}", "false")

    # Extract optimization parameters from the request (if there are any they are in parantheses)
    opt_parameters = np.empty(np.size(request.json["optimizationMethod"]), dtype=object)
    opt_methods = request.json["optimizationMethod"].copy()

    for i in range(np.size(request.json["optimizationMethod"])):
        if "(" in request.json["optimizationMethod"][i]:
            parameters = (
                request.json["optimizationMethod"][i].split("(")[1][:-1].split(",")
            )
            opt_parameters_dict = {}
            for param in parameters:
                key, value = param.split(":")
                opt_parameters_dict[key.strip()] = (
                    float(value) if "." in value else int(value)
                )
            opt_parameters[i] = opt_parameters_dict
            opt_methods[i] = request.json["optimizationMethod"][i].split("(")[0]

    # Enqueue the task with RQ
    job = q.enqueue(
        optimization_function,
        optimization_method=opt_methods,
        latest_design=latest_design,  # This needs to be serializable or reconstructed within the task
        redis_key=job_id,
        additional_opt_parameters=opt_parameters,
        job_timeout="1d",  # Increase the timeout to be able to run long simulations
        # job_id=str(job_id),
    )

    # Store the mapping in Redis or another storage mechanism
    redis_conn.set(f"custom_job_id:{str(job_id)}", job.get_id())

    # Return a success message with the job ID
    return (
        jsonify(
            {
                # "message": "Optimization started successfully",
                # "job_id": str(job_id),
                # "rq_job_id": job.get_id(),
                "merit": initial_merit,
                "step": 0,
            }
        ),
        200,
    )


@optimize_bp.route("/stop_optimization", methods=["POST"])
def stop_optimization():
    job_id = session.get("job_id")
    if not job_id:
        return jsonify({"error": "Job ID is required"}), 400

    # Retrieve the RQ job ID using the custom job ID
    rq_job_id = redis_conn.get(f"custom_job_id:{job_id}")
    if not rq_job_id:
        return jsonify({"error": "Job not found"}), 404

    try:
        job = rq.job.Job.fetch(rq_job_id.decode("utf-8"), connection=redis_conn)
        # If the job is already started, this will raise an exception
        job.cancel()
        # Set a cancellation flag in the redis database that will be retrieved by the actual worker
        redis_conn.set(f"job_cancel_flag:{job_id}", "true")

        return jsonify({"message": "Job cancelled successfully"}), 200
    except (
        rq.exceptions.NoSuchJobError,
        rq.exceptions.InvalidJobOperationError,
        AttributeError,
    ) as e:
        return jsonify({"error": str(e)}), 400


# Now I need a third function that is periodically called via AJAX to update the
# merit graph and the filter representation
@optimize_bp.route("/update_optimization", methods=["GET"])
def update_optimization():
    job_id = session.get("job_id")
    if not job_id:
        # If there's no job_id in the session, return an error
        return jsonify({"error": "Job ID is required"}), 400

    # Check the results from the redis database
    current_best = redis_conn.get(f"current_best:{job_id}")
    current = redis_conn.get(f"current:{job_id}")

    if current_best:
        # Convert the current_best from bytes to a Python dict and return as JSON
        result_data = json.loads(current_best.decode("utf-8"))
        current_merit = json.loads(current.decode("utf-8"))

        # Extract the filter representation from the result data
        (
            num_boxes,
            colors,
            heights,
            number_unique_materials,
            unique_materials,
            unique_colors,
            incoherent,
        ) = extract_filter_design(
            result_data["current_structure"]["structure_thicknesses"],
            result_data["current_structure"]["structure_materials"],
            result_data["current_structure"]["incoherent"],
        )

        # Compile the response object
        response_dict = {
            "best_merit": result_data["merit"],
            "best_step": result_data["step"],
            "current_merit": current_merit["merit"],
            "current_step": current_merit["step"],
            "num_boxes": num_boxes,
            "colors": colors,
            "heights": heights,
            "number_unique_materials": number_unique_materials,
            "unique_materials": unique_materials,
            "unique_colors": unique_colors,
            "incoherent": incoherent,
        }

        # Check if job is already done
        rq_job_id = redis_conn.get(f"custom_job_id:{job_id}")

        if not rq_job_id:
            return jsonify({"error": "Job not found"}), 400

        job = rq.job.Job.fetch(rq_job_id.decode("utf-8"), connection=redis_conn)

        if job.is_finished:
            response_dict["finished"] = True

        # Update database state
        optimization = select_latest_optimization(job_id)
        optimization.current_json = json.dumps(result_data["current_structure"])
        optimization.current_merit = result_data["merit"]
        optimization.steps = result_data["step"]
        db.session.commit()

        return jsonify(response_dict), 200
    else:
        # If no data is available, return a message indicating so
        return jsonify({"message": "No data available"}), 404
