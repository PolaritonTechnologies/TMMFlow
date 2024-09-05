from sqlalchemy import create_engine, Column, Integer, String, Boolean, DateTime
from sqlalchemy.orm import sessionmaker, declarative_base
import pandas as pd
from datetime import datetime

# Database connection setup
DATABASE_URL = (
    "sqlite:///../instance/database.sqlite"  # Replace with your actual database URL
)
engine = create_engine(DATABASE_URL)
Session = sessionmaker(bind=engine)
session = Session()

# Define the base class for declarative models
Base = declarative_base()


class Material(Base):
    __tablename__ = "materials"
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, nullable=False, unique=True)
    creation_time = Column(DateTime, nullable=False, default=datetime.now)
    username = Column(String, nullable=False)
    team = Column(String, nullable=False)
    data = Column(String, nullable=True)


# Create the table if it doesn't exist
Base.metadata.create_all(engine)

# User details
name = "ZrO2"
username = "julian"
# team = "University of Cologne, HCNB"
team = "default"
data = pd.read_csv("../materials/ZrO2.csv", sep="\t", skiprows=1).to_json()

# Create a new material instance
new_material = Material(
    username=username,
    name=name,
    team=team,
    data=data,
)

# Add the new material to the session and commit
session.add(new_material)
session.commit()

print("New material added successfully.")

# Close the session
session.close()
