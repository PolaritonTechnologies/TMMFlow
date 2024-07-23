# Use an official Python runtime as a parent image
FROM python:3.10.12

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt
RUN pip install rq
# Install SQLite3
RUN apt-get update && apt-get install -y sqlite3 && rm -rf /var/lib/apt/lists/*

# Make port 5000 available to the world outside this container
EXPOSE 5000

# Define environment variable pointing to the Flask app in /gui/
ENV FLASK_APP=gui
ENV PYTHONPATH="${PYTHONPATH}:/app/src"
ENV FLASK_ENV=development
ENV FLASK_DEBUG=1

# Run app.py when the container launches
CMD ["flask", "run", "--host=0.0.0.0"]