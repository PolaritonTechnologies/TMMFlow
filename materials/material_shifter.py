import pandas as pd

# Read the CSV file
df = pd.read_csv('C545T.csv', skiprows=1, sep = "\t")

# Subtract 45 from the 'nm' column
df['nm'] = df['nm'] - 45

# Write the data back to the CSV file
df.to_csv('C545T_mod.csv', index=False, sep = "\t")