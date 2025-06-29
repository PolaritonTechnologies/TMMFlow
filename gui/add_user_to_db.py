import bcrypt
from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    Boolean,
    Text,
    Float,
    DateTime,
)
from sqlalchemy.orm import sessionmaker, declarative_base
import json
from datetime import datetime


# Function to hash a password
def hash_password(password):
    salt = bcrypt.gensalt()
    hashed_password = bcrypt.hashpw(password.encode("utf-8"), salt)
    return hashed_password


# Database connection setup
DATABASE_URL = (
    "sqlite:///../instance/database.sqlite"  # Replace with your actual database URL
)
engine = create_engine(DATABASE_URL)
Session = sessionmaker(bind=engine)
session = Session()

# Define the base class for declarative models
Base = declarative_base()


# Define the User model
class User(Base):
    __tablename__ = "users"
    id = Column(Integer, primary_key=True)
    username = Column(String(80), unique=True, nullable=False)
    pw = Column(String(128))
    email = Column(String(80), unique=True, nullable=False)
    team = Column(String(80), nullable=False)
    active = Column(Boolean)

    # Function to verify a password
    def verify_password(self, provided_password):
        return bcrypt.checkpw(provided_password.encode("utf-8"), self.pw)


class Job(Base):
    __tablename__ = "jobs"
    opt_id = Column(Integer, primary_key=True, autoincrement=True)
    job_id = Column(Integer, nullable=True)
    time_stamp = Column(DateTime, nullable=False, default=datetime.now())
    username = Column(String, nullable=False)
    filter_name = Column(String, nullable=False)
    optimization_method = Column(Text, nullable=True)
    current_json = Column(Text, nullable=True)
    current_data = Column(Text, nullable=True)
    steps = Column(Integer, nullable=True)
    initial_merit = Column(Float, nullable=True)
    current_merit = Column(Float, nullable=True)
    description = Column(Text, nullable=True)
    active = Column(Boolean, default=True)


# Create the table if it doesn't exist
Base.metadata.create_all(engine)

# User details
username = "florian"
email = "fleroux@uni-koeln.de"
team = "University of Cologne, HCNB"
password = "TMMFlow"
hashed_password = hash_password(password).decode("utf-8")

# Create a new user instance
new_user = User(
    username=username,
    pw=hashed_password,
    email=email,
    team=team,
    active=True,  # Assuming active is a boolean represented as True or False
)

# Add dummy entry for the new user to the jobs table
# Read in demo_test.json from file
# Open file and read in the JSON
with open("../examples/bandpass_flat_top_600nm.json", "r") as f:
    current_json = json.load(f)

# if the user has no jobs yet, set job_id to 0
if session.query(Job).count() == 0:
    job_id = 0
else:
    latest_job = session.query(Job).order_by(Job.job_id.desc()).first()
    job_id = latest_job.job_id + 1

new_job = Job(
    username=username,
    job_id=job_id,
    filter_name="demo_test",
    current_json=json.dumps(current_json),
)

# Add the new user to the session and commit
session.add(new_user)
session.add(new_job)
session.commit()

print("New user added successfully.")

# Close the session
session.close()
