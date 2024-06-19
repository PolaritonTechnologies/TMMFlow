# The database handled in sqlite by flask need to be modified to add users outside of the webapp running
# open the database.db file with sqlite3 and add the user
# Use python to add users to the database

import sqlite3
import sys

# connect to the database located in instance/database.db
conn = sqlite3.connect("instance/database.db")
c = conn.cursor()

# dump the user table
c.execute("DROP TABLE IF EXISTS users")

# create a table to store users
c.execute(
    "CREATE TABLE IF NOT EXISTS users (id INTEGER PRIMARY KEY AUTOINCREMENT, username TEXT, password TEXT)"
)

usernames = ["Florian", "Julian"]
password = "OugaPw"

for username in usernames:

    # check if the username exist already and then add the user
    c.execute("SELECT * FROM users WHERE username=?", (username,))
    if c.fetchone() is None:
        c.execute(
            "INSERT INTO users (username, password) VALUES (?, ?)", (username, password)
        )
    else:
        print("User already exists")

    # commit the changes
    conn.commit()

# read the current list of users
c.execute("SELECT * FROM users")

# print the list of users but hide the passwords
for user in c.fetchall():
    print(f"User: {user[1]}")

# create a table called jobs which stores the username of the user who ran the job,
# the initial .json file for the job, a list of optimisation performed on the job,
# a .json file for the current job, the number of steps performed on the job initialised at 0,
# the initial merit and the current merit

# delete the job table
# c.execute("DROP TABLE IF EXISTS jobs")

c.execute(
    "CREATE TABLE IF NOT EXISTS jobs (id INTEGER PRIMARY KEY AUTOINCREMENT, username TEXT, initial_json TEXT, optimisations TEXT, current_json TEXT, steps INTEGER, initial_merit REAL, current_merit REAL)"
)

# commit the changes
conn.commit()

# print the job table
c.execute("SELECT * FROM jobs")

print("Jobs: ")
for job in c.fetchall():
    print(f"Job: {job}")

# print steps from the highest jobid
c.execute("SELECT steps FROM jobs ORDER BY id DESC LIMIT 1")
print(f"Steps: {c.fetchone()}")

# close the connection
conn.close()
