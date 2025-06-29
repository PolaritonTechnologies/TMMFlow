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
from sqlalchemy import text


def optimize_db():
    """Configure SQLite for better performance"""
    db.session.execute(text("PRAGMA journal_mode=WAL"))  # Write-Ahead Logging
    db.session.execute(text("PRAGMA synchronous=NORMAL"))  # Faster synchronization
    db.session.execute(text("PRAGMA cache_size=-64000"))  # 64MB cache
    db.session.execute(text("PRAGMA temp_store=MEMORY"))  # Store temp tables in memory
    db.session.commit()


def init_db():
    # If the database exists, you might want to drop and recreate tables or do nothing
    db.drop_all()
    db.create_all()
    optimize_db()


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
    active = db.Column(db.Boolean, default=True)

    # Function to verify a password
    def verify_password(self, provided_password):
        return bcrypt.checkpw(
            provided_password.encode("utf-8"), self.pw.encode("utf-8")
        )


class Job(db.Model):
    __tablename__ = "jobs"
    opt_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    job_id = db.Column(db.Integer, nullable=True)
    time_stamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    username = db.Column(db.String, nullable=False)
    filter_name = db.Column(db.String, nullable=False)
    optimization_method = db.Column(db.Text, nullable=True)
    current_json = db.Column(db.Text, nullable=True)
    current_data = db.Column(db.Text, nullable=True)
    steps = db.Column(db.Integer, nullable=True)
    initial_merit = db.Column(db.Float, nullable=True)
    current_merit = db.Column(db.Float, nullable=True)
    description = db.Column(db.Text, nullable=True)
    active = db.Column(db.Boolean, default=True)

    # Add indexes
    __table_args__ = (
        db.Index("idx_job_id_active", "job_id", "active"),
        db.Index("idx_username_active", "username", "active"),
        db.Index("idx_time_stamp", "time_stamp"),
        db.Index("idx_filter_name", "filter_name"),
        db.Index("idx_job_timestamp_active", "job_id", "time_stamp", "active"),
        db.Index("idx_filter_timestamp_active", "filter_name", "time_stamp", "active"),
    )


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
    creation_time = db.Column(db.DateTime, nullable=False, default=datetime.now)
    username = db.Column(db.String, nullable=False)
    team = db.Column(db.String, nullable=False)
    data = db.Column(db.Text, nullable=True)
    active = db.Column(db.Boolean, default=True)

    # Add indexes
    __table_args__ = (
        db.Index("idx_material_team_active", "team", "active"),
        db.Index("idx_material_name_active", "name", "active"),
    )


def get_latest_job_id():
    return db.session.query(func.max(Job.job_id)).scalar() or 0


def select_latest_optimization(job_id=None, username=None):
    # Use load_only() with actual mapped attributes
    query = Job.query.options(db.load_only(Job.current_json, Job.time_stamp))

    if username:
        return (
            query.filter(Job.username == username, Job.active == True)
            .order_by(Job.time_stamp.desc())
            .first()
        )
    elif job_id != None:
        return (
            query.filter(Job.job_id == job_id, Job.active == True)
            .order_by(Job.time_stamp.desc())
            .first()
        )
    raise Exception("Error loading filter: job_id and username none")


def select_initial_optimization(job_id=None):
    # Load the latest json from the database using the job_id identifier
    if job_id is None:
        raise Exception("Error loading filter: job_id is None")

    # Use indexed columns for filtering
    latest_optimization = (
        Job.query.filter(Job.job_id == job_id, Job.active == True)
        .order_by(Job.time_stamp.asc())
        .first()
    )

    return latest_optimization


def select_optimization_by_datetime(date_time, job_id=None):
    # Truncate the microseconds from the datetime object
    date_time = date_time.replace(microsecond=0)

    # Load the latest json from the database using the job_id identifier and the date_time string that first has to be converted into a datetime object
    filter_status = (
        Job.query.filter((Job.job_id == job_id) & (Job.active == True))
        .filter(
            func.strftime("%Y-%m-%d %H:%M:%S", Job.time_stamp)
            == date_time.strftime("%Y-%m-%d %H:%M:%S")
        )
        .first()
    )

    return filter_status


# Function to delete a specific job from the DB
def deactivate_optimization_from_db(date_time, job_id):
    try:
        filter_status_to_delete = (
            Job.query.options(
                db.load_only(Job.opt_id)
            )  # Only load primary key to deactivate
            .filter(
                Job.job_id == job_id,
                Job.active == True,
                func.strftime("%Y-%m-%d %H:%M:%S", Job.time_stamp)
                == date_time.strftime("%Y-%m-%d %H:%M:%S"),
            )
            .first()
        )

        # Use count() instead of .all() for performance
        number_of_elements = Job.query.filter(
            Job.job_id == job_id, Job.active == True
        ).count()

        if number_of_elements == 1:
            print(
                "Error inactivating filter: Last filter element cannot be deactivated"
            )
            return False

        if filter_status_to_delete:
            filter_status_to_delete.active = False
            db.session.commit()
            return True
        else:
            return False
    except Exception as e:
        db.session.rollback()
        print(f"Error inactivating filter status: {e}")
        return False


def select_job_by_datetime_and_name(date_time, name):
    # Truncate the microseconds from the datetime object
    date_time = date_time.replace(microsecond=0)

    # Use indexed columns for filtering
    filter_status = Job.query.filter(
        Job.filter_name == name,
        Job.active == True,
        func.strftime("%Y-%m-%d %H:%M:%S", Job.time_stamp)
        == date_time.strftime("%Y-%m-%d %H:%M:%S"),
    ).first()

    return filter_status


def deactivate_project_from_db(filter_name, filter_date):
    try:
        with db.session.begin():
            filter_status = select_job_by_datetime_and_name(filter_date, filter_name)
            if not filter_status:
                return False

            db.session.query(Job).filter(
                Job.job_id == filter_status.job_id, Job.active == True
            ).update({"active": False})

        return True
    except Exception as e:
        db.session.rollback()
        print(f"Error deactivating project: {e}")
        return False


def load_latest_filter(job_id=None, username=None):
    # Load the latest filter from the database using the job_id identifier
    if job_id is None and username is not None:
        latest_optimization = select_latest_optimization(username=username)
        if latest_optimization and latest_optimization.current_json:
            latest_json = json.loads(latest_optimization.current_json)
        else:
            raise Exception("No optimization found for the given username.")
    elif job_id is not None:
        latest_optimization = select_latest_optimization(job_id=job_id)
        if latest_optimization and latest_optimization.current_json:
            latest_json = json.loads(latest_optimization.current_json)
        else:
            raise Exception("No optimization found for the given job_id.")
    else:
        raise Exception("Error loading filter: job_id and username are None")

    # Now create the FilterStack object using the initial json and modify so
    # that it matches the latest json (needed to comply with the FilterStack
    # logic)
    my_filter = FilterStack(my_filter_dict=latest_json)

    return my_filter


def get_all_filter_versions(job_id=None):
    return (
        Job.query.options(db.load_only(Job.current_json, Job.time_stamp))
        .filter(Job.job_id == job_id, Job.active == True)
        .order_by(Job.time_stamp.desc())
        .all()
    )


def get_all_user_projects(user=None):
    # Optimize with single query using window functions and indexed columns
    stmt = (
        db.session.query(
            Job.job_id,
            Job.filter_name,
            Job.time_stamp,
            func.row_number()
            .over(partition_by=Job.job_id, order_by=Job.time_stamp.desc())
            .label("rn"),
        )
        .filter(Job.username == user, Job.active == True)
        .subquery()
    )

    results = (
        db.session.query(stmt.c.job_id, stmt.c.filter_name, stmt.c.time_stamp)
        .filter(stmt.c.rn == 1)
        .all()
    )

    return (
        [r.job_id for r in results],
        [r.filter_name for r in results],
        [r.time_stamp for r in results],
    )


# Function to hash a password
def change_password_database(username, password):
    # Select user using indexed columns
    user = User.query.filter(User.username == username, User.active == True).first()

    if not user:
        raise Exception("User not found or inactive.")

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
        ((Material.team == team) | (Material.team == "default"))
        & (Material.active == True)
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
        Material.name == material,
        Material.active == True,
        (Material.team == team) | (Material.team == "default"),
    ).first()

    if not material_data or not material_data.data:
        raise Exception("Material not found or no data available.")

    # Convert material_data from json to pandas dataframe
    material_data_df = pd.DataFrame(json.loads(material_data.data))

    return material_data_df


# Function to deactivate a specific material from the DB
def deactivate_material_from_db(material, team):
    try:
        # Find the material to deactivate using indexed columns
        material_to_delete = Material.query.filter(
            Material.name == material, Material.team == team, Material.active == True
        ).first()

        if material_to_delete:
            # Deactivate the material
            material_to_delete.active = False
            # Commit the changes
            db.session.commit()
            return True
        else:
            return False
    except Exception as e:
        db.session.rollback()
        print(f"Error deactivating material: {e}")
        return False


def add_material_to_db(material, team, material_class, material_data, username):
    # Create a new material object
    new_material = Material(
        name=material,
        material_class=material_class,
        username=username,
        team=team,
        data=material_data,
        active=True,
    )

    # Add the new material to the database
    db.session.add(new_material)

    # Commit the changes
    db.session.commit()
    return True
