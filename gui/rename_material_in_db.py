import sqlite3
import pandas as pd
import json

# Connect to the SQLite database
conn = sqlite3.connect("./instance/database.sqlite")
cursor = conn.cursor()

# Retrieve all rows from the materials table
cursor.execute("SELECT rowid, data FROM materials")
rows = cursor.fetchall()


# Function to rename columns based on the number of columns
def rename_columns(df):
    empty_dict = {}
    if len(df.columns) == 3:
        df.columns = ["wavelength", "n", "k"]
    elif len(df.columns) == 5:
        df.columns = ["wavelength", "n_ord", "k_ord", "n_exord", "k_exord"]
    elif len(df.columns) == 7:
        df.columns = ["wavelength", "n_x", "k_x", "n_y", "k_y", "n_z", "k_z"]

    for col in df.columns:
        empty_dict[col] = df[col].tolist()

    return empty_dict


# Process each row
for row in rows:
    rowid = row[0]
    json_data = row[1]

    # Convert JSON data to pandas DataFrame
    df = pd.read_json(json_data)

    # Rename columns
    empty_dict = rename_columns(df)

    # Convert the reformatted JSON back to a string
    new_json_data = json.dumps(empty_dict)

    # pd.DataFrame.from_dict(json.loads(new_json_data))

    # Update the database with the new JSON data
    cursor.execute(
        "UPDATE materials SET data = ? WHERE rowid = ?", (new_json_data, rowid)
    )

# Commit the changes and close the connection
conn.commit()
conn.close()
