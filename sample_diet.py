import pandas as pd
import numpy as np


def sample_diet(diet_csv_path, variability=0.1):
    """
    Randomly samples metabolite amounts from a dietary composition file,
    assuming fixed nominal amounts with optional variability.

    Args:
        diet_csv_path (str): Path to a CSV file with columns:
            'Metabolite ID', 'Amount (mmol)'.
        variability (float): Percentage variability to simulate dietary fluctuation.
                             Default is 0.1 (i.e., ±10%).

    Returns:
        dict: A dictionary where keys are metabolite IDs and values are
              randomly perturbed amounts in mmol.
    """

    df = pd.read_csv(diet_csv_path)

    df["Sampled Amount (mmol)"] = df["Amount (mmol)"] * np.random.uniform(
        1 - variability, 1 + variability, size=len(df)
    )

    sampled_amounts = dict(zip(df["Metabolite ID"], df["Sampled Amount (mmol)"]))

    return sampled_amounts


def sample_gases(T=310, R=0.08206, P=1):
    """
    Samples gas volumes from normal distributions and converts to mmol.

    Args:
        T (float): Temperature in Kelvin (default: 310 K).
        R (float): Ideal gas constant in L·atm/(mol·K) (default: 0.08206).
        P (float): Total pressure in atm (default: 1 atm).

    Returns:
        dict: Mapping from gas metabolite ID to amount in mmol.
    """

    gas_volume_stats = {
        "o2[e]": (0.58, 0.43),
        "co2[e]": (9.7, 2.4),
        "n2[e]": (64, 52),
        "h2[e]": (14, 9.9),
        "ch4[e]": (5.6, 7.6)
    }

    # Sample individual gas volumes in mL
    sampled_volumes_ml = {
        gas: max(0, np.random.normal(mean, std))  # Clip to avoid negative volumes
        for gas, (mean, std) in gas_volume_stats.items()
    }

    # Convert total volume to litres
    total_volume_L = sum(sampled_volumes_ml.values()) / 1000

    # Compute v/v % for each gas
    vv_fractions = {
        gas: vol_ml / sum(sampled_volumes_ml.values())
        for gas, vol_ml in sampled_volumes_ml.items()
    }

    # Apply ideal gas law to get mmol
    gas_mmol = {
        gas: 1000 * frac * P * total_volume_L / (R * T)
        for gas, frac in vv_fractions.items()
    }

    return gas_mmol
