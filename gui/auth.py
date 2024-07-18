from flask import (
    current_app,
    Blueprint,
    flash,
    redirect,
    render_template,
    request,
    session,
    url_for,
    jsonify,
)
from flask_login import (
    LoginManager,
    login_required,
    login_user,
    logout_user,
    current_user,
)
from . import db
from .database import User

auth_bp = Blueprint("auth", __name__, url_prefix="/auth")

login_manager = LoginManager()
# Adjusted to include the blueprint name
login_manager.login_view = "auth.login"


@login_manager.user_loader
def load_user(user_id):
    return User.query.get(user_id)


# The clear_session function might need to be adjusted or moved depending on your app structure
@auth_bp.before_app_request
def clear_session():
    db.create_all()
    # Check if clear_session is in the list before attempting to remove it
    if clear_session in current_app.before_request_funcs.get(None, []):
        current_app.before_request_funcs[None].remove(clear_session)
    session.clear()


@auth_bp.context_processor
def inject_user():
    return dict(logged_in_user=current_user)


@auth_bp.route("/logout")
@login_required
def logout():
    logout_user()
    session.clear()
    return redirect(url_for("auth.login"))  # Adjusted to include the blueprint name


@auth_bp.route("/login", methods=("GET", "POST"))
def login():
    if request.method == "POST":
        username = request.form["username"]
        password = request.form["password"]
        user = User.query.filter_by(username=username).first()

        if user is not None and user.check_password(password):
            session["user_id"] = username
            session["job_id"] = None
            session["filter_name"] = None
            login_user(user)
            flash("Logged in successfully.")
            return redirect(url_for("stack_bp.stack"))
        else:
            flash("Invalid username or password.")

    return render_template("login.html")


@auth_bp.route("/handle_send_username", methods=["POST"])
def handle_send_username():
    username = request.json["username"]
    session["user_id"] = username

    return jsonify({"message": f"Username {username} received"})
