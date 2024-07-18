# Flask related imports
from flask import current_app
from flask.cli import with_appcontext
import click
from flask_login import UserMixin

# Standard pyton functions
from datetime import datetime
import json
import os

# Import the SQLAlchemy instance
from . import db

# Import the FilterStack class
from FilterStack import FilterStack


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
    __tablename__ = "users"
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    password = db.Column(db.String(128))

    def check_password(self, password):
        return password == self.password


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


def get_latest_job_id():
    # Get the latest job_id from the database
    latest_job_id = Job.query.order_by(Job.job_id.desc()).first()

    if latest_job_id == None:
        return 0
    else:
        return int(latest_job_id.job_id)


def select_latest_optimization(job_id=None):
    # Load the latest json from the database using the job_id identifier
    if job_id is None:
        raise Exception("Error loading filter: job_id is None")

    # Now execute an SQL query to look for all optimization entries in
    # optimizations table with the job_id and pick the one with the latest
    # time_stamp value
    latest_optimization = (
        Job.query.filter_by(job_id=job_id).order_by(Job.time_stamp.desc()).first()
    )

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


def load_latest_filter(job_id=None):
    # Load the latest filter from the database using the job_id identifier
    if job_id is None:
        raise Exception("Error loading filter: job_id is None")

    # Now execute an SQL query to look for all optimization entries in
    # optimizations table with the job_id and pick the one with the latest
    # time_stamp value
    latest_json = json.loads(select_latest_optimization(job_id).current_json)

    # Now create the FilterStack object using the initial json and modify so
    # that it matches the latest json (needed to comply with the FilterStack
    # logic)
    my_filter = FilterStack(my_filter_dict=latest_json)

    return my_filter
