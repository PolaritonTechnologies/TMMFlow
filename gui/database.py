# Flask related imports
from flask import current_app
from flask.cli import with_appcontext
import click
from flask_login import UserMixin

# Standard python functions
from datetime import datetime
import json
import pandas as pd

# Hash passwords for security
import bcrypt

# Import the SQLAlchemy instance
from . import db

# Import the FilterStack class
from FilterStack import FilterStack
from sqlalchemy import func, distinct


def init_db():
    # If the database exists, you might want to drop and recreate tables or do nothing
    db.drop_all()
    db.create_all()


@click.command("init-db")
@with_appcontext
def init_db_command():
    """Clear existing data and create new tables."""
    init_db()
    click.echo("Initialized the database.")


# Define the database models using the module-level db instance
class User(db.Model, UserMixin):
    """
    CREATE TABLE users (
    id INTEGER PRIMARY KEY,
    username TEXT,
    pw TEXT,
    email TEXT,
    team TEXT,
    active BOOLEAN
    );
    """

    __tablename__ = "users"
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    pw = db.Column(db.String(128))
    email = db.Column(db.String(80), unique=True, nullable=False)
    team = db.Column(db.String(80), nullable=False)
    active = db.Column(db.Boolean)

    # Function to verify a password
    def verify_password(self, provided_password):
        return bcrypt.checkpw(
            provided_password.encode("utf-8"), self.pw.encode("utf-8")
        )


class Job(db.Model):
    __tablename__ = "jobs"
    opt_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    job_id = db.Column(db.Integer, nullable=True)
    time_stamp = db.Column(db.DateTime, nullable=False, default=datetime.now())
    username = db.Column(db.String, nullable=False)
    filter_name = db.Column(db.String, nullable=False)
    optimization_method = db.Column(db.Text, nullable=True)
    current_json = db.Column(db.Text, nullable=True)
    current_data = db.Column(db.Text, nullable=True)
    steps = db.Column(db.Integer, nullable=True)
    initial_merit = db.Column(db.Float, nullable=True)
    current_merit = db.Column(db.Float, nullable=True)


class Material(db.Model):
    """
    CREATE TABLE materials (
    id INTEGER PRIMARY KEY,
    name TEXT,
    creation_time DATETIME DEFAULT CURRENT_TIMESTAMP,
    username TEXT,
    team TEXT,
    data TEXT
    );
    """

    __tablename__ = "materials"
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String, nullable=False, unique=True)
    material_class = db.Column(db.String, nullable=False)
    creation_time = db.Column(db.DateTime, nullable=False, default=datetime.now())
    username = db.Column(db.String, nullable=False)
    team = db.Column(db.String, nullable=False)
    data = db.Column(db.Text, nullable=True)


def get_latest_job_id():
    # Get the latest job_id from the database
    latest_job_id = Job.query.order_by(Job.job_id.desc()).first()

    if latest_job_id == None:
        return 0
    else:
        return int(latest_job_id.job_id)


def select_latest_optimization(job_id=None, username=None):
    # Load the latest json from the database using the job_id identifier
    if username is not None and job_id is None:
        # Get the latest Job for the given user overall
        latest_optimization = (
            Job.query.filter_by(username=username)
            .order_by(Job.time_stamp.desc())
            .first()
        )
    elif job_id is not None:
        # Now execute an SQL query to look for all optimization entries in
        # optimizations table with the job_id and pick the one with the latest
        # time_stamp value
        latest_optimization = (
            Job.query.filter_by(job_id=job_id).order_by(Job.time_stamp.desc()).first()
        )
    else:
        raise Exception("Error loading filter: job_id and username none")

    return latest_optimization


def select_initial_optimization(job_id=None):
    # Load the latest json from the database using the job_id identifier
    if job_id is None:
        raise Exception("Error loading filter: job_id is None")

    # Now execute an SQL query to look for all optimization entries in
    # optimizations table with the job_id and pick the one with the latest
    # time_stamp value
    latest_optimization = (
        Job.query.filter_by(job_id=job_id).order_by(Job.time_stamp.asc()).first()
    )

    return latest_optimization


def select_optimization_by_datetime(date_time, job_id=None):
    # Truncate the microseconds from the datetime object
    date_time = date_time.replace(microsecond=0)

    # Load the latest json from the database using the job_id identifier and the date_time string that first has to be converted into a datetime object
    filter_status = (
        Job.query.filter_by(job_id=job_id)
        .filter(
            func.strftime("%Y-%m-%d %H:%M:%S", Job.time_stamp)
            == date_time.strftime("%Y-%m-%d %H:%M:%S")
        )
        .first()
    )

    return filter_status


def select_job_by_datetime_and_name(date_time, name):
    # Truncate the microseconds from the datetime object
    date_time = date_time.replace(microsecond=0)

    # Load the latest json from the database using the job_id identifier and the date_time string that first has to be converted into a datetime object
    filter_status = (
        Job.query.filter_by(filter_name=name)
        .filter(
            func.strftime("%Y-%m-%d %H:%M:%S", Job.time_stamp)
            == date_time.strftime("%Y-%m-%d %H:%M:%S")
        )
        .first()
    )

    return filter_status


def load_latest_filter(job_id=None, username=None):
    # Load the latest filter from the database using the job_id identifier
    if job_id is None and username is not None:
        latest_json = json.loads(
            select_latest_optimization(username=username).current_json
        )
    elif job_id is not None:
        latest_json = json.loads(select_latest_optimization(job_id).current_json)
    else:
        raise Exception("Error loading filter: job_id and username are None")

    # Now create the FilterStack object using the initial json and modify so
    # that it matches the latest json (needed to comply with the FilterStack
    # logic)
    my_filter = FilterStack(my_filter_dict=latest_json)

    return my_filter


def get_all_filter_versions(job_id=None):
    # Load all filter versions from the database
    all_filter_versions = (
        Job.query.filter_by(job_id=job_id).order_by(Job.time_stamp.desc()).all()
    )

    return all_filter_versions


def get_all_user_projects(user=None):
    # Load all unique job IDs from the database for the given user
    unique_job_ids = (
        Job.query.with_entities(distinct(Job.job_id)).filter_by(username=user).all()
    )

    # Extract job IDs from the result
    unique_job_ids = [job_id[0] for job_id in unique_job_ids]

    # Load filter names and timestamps for the unique job IDs
    filter_names = []
    time_stamps = []
    for job_id in unique_job_ids:
        filter_name, time_stamp = (
            Job.query.with_entities(Job.filter_name, Job.time_stamp)
            .filter_by(job_id=job_id)
            .order_by(Job.time_stamp.desc())
            .first()
        )
        if filter_name and time_stamp:
            filter_names.append(filter_name)
            time_stamps.append(time_stamp)

    return unique_job_ids, filter_names, time_stamps


def get_all_user_projects(user=None):
    # Load all unique job IDs from the database for the given user
    unique_job_ids = (
        Job.query.with_entities(distinct(Job.job_id)).filter_by(username=user).all()
    )

    # Extract job IDs from the result
    unique_job_ids = [job_id[0] for job_id in unique_job_ids]

    # Load filter names and timestamps for the unique job IDs
    filter_names = []
    time_stamps = []
    for job_id in unique_job_ids:
        filter_name, time_stamp = (
            Job.query.with_entities(Job.filter_name, Job.time_stamp)
            .filter_by(job_id=job_id)
            .order_by(Job.time_stamp.desc())
            .first()
        )
        if filter_name and time_stamp:
            filter_names.append(filter_name)
            time_stamps.append(time_stamp)

    return unique_job_ids, filter_names, time_stamps


# Function to hash a password
def change_password_database(username, password):
    # Select user
    user = User.query.filter_by(username=username).first()

    # Generate a salt
    salt = bcrypt.gensalt()
    # Hash the password with the salt
    hashed_password = bcrypt.hashpw(password.encode("utf-8"), salt)

    # Update the user's password
    user.pw = hashed_password.decode("utf-8")

    # Write back to database
    db.session.commit()


# Function to get all available materials from DB
def get_available_materials(team):
    # Get all available materials for the given team and the default team "default"
    available_materials_db_entries = Material.query.filter(
        (Material.team == team) | (Material.team == "default")
    ).all()

    # Extract the names of the materials
    available_materials = [material.name for material in available_materials_db_entries]
    material_classes = [
        material.material_class for material in available_materials_db_entries
    ]

    return available_materials, material_classes


# Function to get a specific material from DB
def get_material_data(material, team):
    # Get the material data for the specific material "material" and the team "team"
    material_data = Material.query.filter(
        (Material.name == material)
        & ((Material.team == team) | (Material.team == "default"))
    ).first()

    # Convert material_data from json to pandas dataframe
    material_data_df = pd.DataFrame(json.loads(material_data.data))

    return material_data_df


# Function to delete a specific material from the DB
def delete_material_from_db(material, team):
    try:
        # Find the material to delete
        material_to_delete = Material.query.filter_by(name=material, team=team).first()

        if material_to_delete:
            # Delete the material
            db.session.delete(material_to_delete)
            # Commit the changes
            db.session.commit()
            return True
        else:
            return False
    except Exception as e:
        db.session.rollback()
        print(f"Error deleting material: {e}")
        return False
