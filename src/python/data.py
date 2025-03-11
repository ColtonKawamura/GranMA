import pandas as pd
import numpy as np


def filter_data(df, *args):
    filtered_df = df.copy()
    
    for i in range(0, len(args), 2):
        value = args[i]
        column_name = args[i+1]
        
        if column_name in filtered_df.columns:
            # Find the row with the closest value
            closest_index = (filtered_df[column_name] - value).abs().idxmin()
            closest_value = filtered_df.loc[closest_index, column_name]
            
            # Filter DataFrame based on the closest value
            filtered_df = filtered_df[filtered_df[column_name] == closest_value]
        else:
            print(f"Warning: Column '{column_name}' not found in DataFrame. Skipping filter.")

    return filtered_df