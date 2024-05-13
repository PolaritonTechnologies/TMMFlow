import matplotlib.pyplot as plt
import os

def allowed_file(filename):
    ALLOWED_EXTENSIONS = {"json", "ofp"}
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def generate_colors(n):
    cmap = plt.cm.get_cmap("viridis", n)  # Get the 'viridis' color map
    colors = [cmap(i) for i in range(cmap.N)]  # Generate colors
    # Convert RGB colors to hex
    hex_colors = [
        "#" + "".join([format(int(c * 255), "02X") for c in color[:3]])
        for color in colors
    ]
    return hex_colors


def get_available_materials(directory):
    material_list = [
        os.path.splitext(f)[0]
        for f in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, f)) and f.endswith(".csv")
    ]
    return material_list

def is_number(string):
    try:
        int(string)
        return True
    except ValueError:
        return False