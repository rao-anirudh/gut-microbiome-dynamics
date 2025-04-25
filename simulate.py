from utilities import *
from sample_diet import sample_diet, sample_gases
from sample_phyla import sample_microbial_library
import warnings
import logging
import time
import os
import pandas as pd
from datetime import datetime

# Suppress warnings during execution
warnings.filterwarnings("ignore")

# Record the simulation start time for unique file naming
sim_time = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')


# Function to record metabolome data into a CSV file
def record_metabolome(t, metabolome, filename):
    """
    Appends the metabolome data at time 't' to the specified file.

    Parameters:
    - t (int): Time point of data collection.
    - metabolome (dict): Dictionary containing metabolite data.
    - filename (str): Path to the CSV file where data will be recorded.
    """
    metabolome = pd.DataFrame({str(t): metabolome})

    # Check if file already exists
    if os.path.exists(filename):
        df = pd.read_csv(filename, index_col=0)
        df = df.join(metabolome, how='outer')
    else:
        df = metabolome

    # Sort the columns by time (ascending)
    df = df[sorted(df.columns, key=int)]
    df.to_csv(filename)


# Function to record microbiome data into a CSV file
def record_microbiome(t, microbiome, filename):
    """
    Appends the microbiome data at time 't' to the specified file.

    Parameters:
    - t (int): Time point of data collection.
    - microbiome (dict): Dictionary containing microbiome species data.
    - filename (str): Path to the CSV file where data will be recorded.
    """
    microbiome = pd.DataFrame({str(t): microbiome})

    # Check if file already exists
    if os.path.exists(filename):
        df = pd.read_csv(filename, index_col=0)
        df = df.join(microbiome, how='outer')
    else:
        df = microbiome

    # Sort the columns by time (ascending)
    df = df[sorted(df.columns, key=int)]

    df.to_csv(filename)


# Function to record growth rate data into a CSV file
def record_growth_rate(t, growth_rate, filename):
    """
    Appends the growth rate data at time 't' to the specified file.

    Parameters:
    - t (int): Time point of data collection.
    - growth_rate (float): Growth rate of the microbiome.
    - filename (str): Path to the CSV file where data will be recorded.
    """
    growth_rate_df = pd.DataFrame({str(t): [growth_rate]})

    # Check if file already exists
    if os.path.exists(filename):
        df = pd.read_csv(filename, index_col=0)
        df = df.join(growth_rate_df, how='outer')
    else:
        df = growth_rate_df

    # Sort the columns by time (ascending)
    df = df[sorted(df.columns, key=int)]
    df.to_csv(filename)


# Main simulation function
def simulate(duration, diet_file, seed=5240):
    """
    Simulates the gut microbiome and metabolome over a specified duration.

    Parameters:
    - duration (int): Total duration of the simulation (in hours).
    - diet_file (str): Path to the diet CSV file to sample diet data.
    - seed (int): Random seed for reproducibility.
    """
    np.random.seed(seed)  # Set random seed for reproducibility

    logging.getLogger("cobra").setLevel(logging.ERROR)  # Suppress cobra library warnings

    # Extract diet name from the diet file for folder naming
    diet_name = diet_file.split("_")[0]
    folder_name = f"{sim_time}_{diet_name}"

    # Create directory to store results
    results_dir = os.path.join("results", folder_name)
    os.makedirs(results_dir, exist_ok=True)

    # Define file paths for saving data
    small_intestine_metabolome_file = os.path.join(results_dir,
                                                   f"{sim_time}_{diet_name}_small_intestine_metabolome.csv")
    small_intestine_microbiome_file = os.path.join(results_dir,
                                                   f"{sim_time}_{diet_name}_small_intestine_microbiome.csv")
    large_intestine_metabolome_file = os.path.join(results_dir,
                                                   f"{sim_time}_{diet_name}_large_intestine_metabolome.csv")
    large_intestine_microbiome_file = os.path.join(results_dir,
                                                   f"{sim_time}_{diet_name}_large_intestine_microbiome.csv")
    small_intestine_growth_file = os.path.join(results_dir, f"{sim_time}_{diet_name}_small_intestine_growth.csv")
    large_intestine_growth_file = os.path.join(results_dir, f"{sim_time}_{diet_name}_large_intestine_growth.csv")

    # Instantiate small and large intestine objects
    small_intestine = SmallIntestine()
    large_intestine = LargeIntestine()

    # Start simulation time
    t = 0

    # Run the simulation for the specified duration
    while t < duration:

        si_growth_rates = dict()
        li_growth_rates = dict()

        # Simulate for small intestine at intervals based on input frequency
        if t % small_intestine.input_frequency == 0:
            print(t)

            # Sample diet and gases and update the small intestine
            sampled_diet = sample_diet(diet_file)
            sampled_gases = sample_gases()
            sampled_diet.update(sampled_gases)
            small_intestine.add_to_metabolome(sampled_diet)

            # Sample microbial library and add to small intestine microbiome
            sampled_microbes = sample_microbial_library("representative_strains.json")
            small_intestine.add_to_microbiome(sampled_microbes)

            # Simulate metabolism and get growth rates for the small intestine
            si_growth_rates = small_intestine.metabolise()

            t += small_intestine.output_frequency  # Update time by the small intestine output frequency

            # Record data for small intestine
            record_microbiome(t, small_intestine.microbiome, small_intestine_microbiome_file)
            record_metabolome(t, small_intestine.metabolome, small_intestine_metabolome_file)
            record_growth_rate(t, small_intestine.growth_rate, small_intestine_growth_file)

        # Simulate transfer from small intestine to large intestine at specific time intervals
        if t % small_intestine.input_frequency == 4:
            print(t)

            small_intestine.transfer(large_intestine, si_growth_rates)

            # Simulate metabolism for the large intestine
            li_growth_rates = large_intestine.metabolise()

            t += large_intestine.output_frequency - 4  # Update time by the large intestine output frequency

            # Record data for large intestine
            record_microbiome(t, large_intestine.microbiome, large_intestine_microbiome_file)
            record_metabolome(t, large_intestine.metabolome, large_intestine_metabolome_file)
            record_growth_rate(t, large_intestine.growth_rate, large_intestine_growth_file)

        # Simulate further transfer and interactions within large intestine at specific intervals
        if t % large_intestine.output_frequency == 0:
            print(t)

            large_intestine.transfer(li_growth_rates)


# Entry point of the simulation
if __name__ == "__main__":
    from multiprocessing import freeze_support

    freeze_support()  # Freeze support for multiprocessing
    start = time.time()  # Record start time for performance measurement
    simulate(24 * 365, "keto_diet.csv")  # Run the simulation for 1 year with keto diet
    stop = time.time()  # Record end time

    # Print the total time taken for the simulation
    print(f"\nTime taken = {(stop - start) / 60} minutes")
