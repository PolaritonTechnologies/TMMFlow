import pandas as pd

material_name = "Ta2O5"
# Read the file
df = pd.read_csv(f"{material_name}.mat", sep=",", header=None)

# Rename the columns
df.columns = ["nm", "n", "k"]

# Add a header row
header = pd.DataFrame(
    [[f"Opt. Const. of {material_name} vs. nm", "", ""], ["nm", "n", "k"]],
    columns=df.columns,
)
df = pd.concat([header, df])

# Write the file
df.to_csv(f"{material_name}.csv", index=False, header=False)
