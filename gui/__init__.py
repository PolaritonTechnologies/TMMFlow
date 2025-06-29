# Flask related imports
from flask import Flask, render_template
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.orm import DeclarativeBase
from redis import Redis
import rq

# import rq_dashboard

# Standard pyton functions
import os


# Define database globally
class Base(DeclarativeBase):
    pass


db = SQLAlchemy(model_class=Base)

# Define the Redis connection and the RQ queue
redis_url = os.environ.get(
    "REDIS_URL", "redis://localhost:6379"
)  # Default to localhost if not set
redis_conn = Redis.from_url(redis_url)
# q = rq.Queue(connection=redis_conn)
# redis_conn = Redis()
# q = rq.Queue(connection=redis_conn)


def create_app(test_config=None):
    """
    For how to do this with SQLAlchemy check out
    https://github.com/pallets-eco/flask-sqlalchemy/blob/main/examples/flaskr/flaskr/__init__.py
    """
    app = Flask(__name__, instance_relative_config=True)

    # some deploy systems set the database url in the environ
    db_url = os.environ.get("DATABASE_URL")

    if db_url is None:
        # default to a sqlite database in the instance folder
        db_url = "sqlite:///database.sqlite"

    app.config.from_mapping(
        # default secret that should be overridden in environ or config
        SECRET_KEY=os.environ.get("SECRET_KEY", "dev"),
        SQLALCHEMY_DATABASE_URI=db_url,
        SQLALCHEMY_POOL_SIZE=10,
        SQLALCHEMY_POOL_TIMEOUT=30,
        SQLALCHEMY_POOL_RECYCLE=1800,
        SQLALCHEMY_MAX_OVERFLOW=20,
        UPLOAD_FOLDER="./instance/temp/",
        MATERIAL_FOLDER="./materials/",
        TEMPLATE_FOLDER="./examples/",
        DEFAULT_FILE="./examples/demo_test.json",
        SESSION_FILES="./gui/session_folder/",
        SESSION_TYPE="filesystem",
        SESSION_COOKIE_SAMESITE="Lax",
        # RQ_DASHBOARD_REDIS_URL="redis://localhost:6379",
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile("config.py", silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder and within it the temp folder exist.
    # os.makedirs fails if it does already exist
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass
    try:
        os.makedirs(app.instance_path + "/temp/")
    except OSError:
        pass

    # initialize Flask-SQLAlchemy and the init-db command
    from .database import init_db_command

    db.init_app(app)
    app.cli.add_command(init_db_command)

    from . import auth, stack, simulate, materials, optimize, settings

    auth.login_manager.init_app(app)

    app.register_blueprint(auth.auth_bp, url_prefix="/")
    app.register_blueprint(stack.stack_bp, url_prefix="/")
    app.register_blueprint(simulate.simulate_bp, url_prefix="/")
    app.register_blueprint(materials.materials_bp, url_prefix="/")
    app.register_blueprint(optimize.optimize_bp, url_prefix="/")
    app.register_blueprint(settings.settings_bp, url_prefix="/")

    # app.add_url_rule("/", endpoint="stack_bp.stack")
    @app.route("/")
    def start_page():
        return render_template("startpage.html")

    # app.config.from_object(rq_dashboard.default_settings)
    # rq_dashboard.web.setup_rq_connection(app)
    # app.register_blueprint(rq_dashboard.blueprint, url_prefix="/rq")

    return app
