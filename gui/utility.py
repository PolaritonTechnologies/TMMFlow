import matplotlib.pyplot as plt

def allowed_file(filename):
    ALLOWED_EXTENSIONS = {"json"}
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
