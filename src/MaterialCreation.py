import pandas as pd
import numpy as np


def constant_index(refractive_index, name):
    """
    generate a constant index of a material.

    Parameters
    ----------
    refractive_index : float
        The refractive index of the material.
    """

    # Create a DataFrame
    df = pd.DataFrame(
        {"nm": range(200, 1601), "n": [1] * (1601 - 200), "k": [0] * (1601 - 200)}
    )

    # Add a header row
    header = pd.DataFrame(
        [[f"Opt. Const. of {name} vs. nm", "", ""], ["nm", "n", "k"]],
        columns=df.columns,
    )
    df = pd.concat([header, df])

    # Write the file
    df.to_csv(f"{name}.csv", index=False, header=False)


def Cauchy_index(A, B, C, name, Ak=0.0, exponent=1.0, edge=4000):
    """
    generate a Cauchy index of a material.

    Parameters
    ----------
    A : float
        The A coefficient of the Cauchy index.
    B : float
        The B coefficient of the Cauchy index.
    C : float
        The C coefficient of the Cauchy index.
    Ak : float
        The Ak coefficient of the Cauchy index.
    exponent : float
        The exponent of the Cauchy index.
    edge : float
        The edge of the Cauchy index.
    """

    # Generate the wavelengths
    wavelengths = np.arange(200, 1601)  # in nm

    # Convert wavelengths to microns
    wavelengths_micron = wavelengths * 0.001

    # Calculate the refractive indices
    n_values = A + B / wavelengths_micron**2 + C / wavelengths_micron**4
    k_values = -Ak * np.exp(
        12400.0 * exponent * ((1.0 / (10000.0 * wavelengths_micron)) - (1.0 / edge))
    )

    # Create a DataFrame
    df = pd.DataFrame({"nm": wavelengths, "n": n_values, "k": k_values})

    # Add a header row
    header = pd.DataFrame(
        [[f"Opt. Cauchy of {name} vs. nm", "", ""], ["nm", "n", "k"]],
        columns=df.columns,
    )
    df = pd.concat([header, df])

    # Write the file
    df.to_csv(f"{name}.csv", index=False, header=False)


# Fused Silica
Cauchy_index(1.4466, 0.0064509, -0.00015588, "FusedSilica")
