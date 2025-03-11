import pandas as pd
import numpy as np

from src.python.data import filter_data
from src.python.plotting import plot_amp

# Load the CSV file
df = pd.read_csv("out/processed/2d_bi_K100_W5.csv")

filtered_df = filter_data(df, 0.1, "pressure", 0.1, "omega", 0.05, "gamma", 1, "seed")

result = plot_amp(filtered_df, plot=True, shear=False, transverse_axis="y")
print("Result:", result)