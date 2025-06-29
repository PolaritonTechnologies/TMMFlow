import sqlite3
import pandas as pd
import json

# Connect to the SQLite database
conn = sqlite3.connect("../instance/database.sqlite")
cursor = conn.cursor()

# Retrieve all rows from the materials table
opt_id = "1249"


# Set active to 0 for the selected row
cursor.execute("UPDATE jobs SET active = ? WHERE opt_id = ?", (0, opt_id))

# Commit the changes and close the connection
conn.commit()
conn.close()
