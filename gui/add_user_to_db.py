import bcrypt
from sqlalchemy import create_engine, Column, Integer, String, Boolean
from sqlalchemy.orm import sessionmaker, declarative_base


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


# Create the table if it doesn't exist
Base.metadata.create_all(engine)

# User details
username = "julian"
email = "julian.butscher@uni-koeln.de"
team = "University of Cologne, HCNB"
password = "OugaPw"
hashed_password = hash_password(password).decode("utf-8")

# Create a new user instance
new_user = User(
    username=username,
    pw=hashed_password,
    email=email,
    team=team,
    active=True,  # Assuming active is a boolean represented as True or False
)

# Add the new user to the session and commit
session.add(new_user)
session.commit()

print("New user added successfully.")

# Close the session
session.close()
