import pandas as pd
import math
import numpy as np

# Load the CSV files into DataFrames
df1 = pd.read_csv('prediction.csv')
df2 = pd.read_csv('output_vdm.csv')
df3 = pd.read_csv('training.csv')

# Specify the common column
common_column = 'TIME'  # Replace with the name of the common column

# Merge the DataFrames based on the common column
merged_df = pd.merge(df1, df2, on=common_column)
# merged_df = pd.concat([df1,df2])
merged_df2 = pd.merge(df3, df2, on=common_column)

# merged_df.loc[np.sqrt(merged_df['VX']**2 + merged_df['VY']**2) <= 0.05, 'AX'] = 0
# merged_df.loc[np.sqrt(merged_df['VX']**2 + merged_df['VY']**2) <= 0.05, 'AY'] = 0

# Save the resulting DataFrame to a new CSV file
merged_df.to_csv('predict_compare.csv', index=False)
merged_df2.to_csv('training_data.csv')




print("Rows with common values have been merged and saved to 'merged_file.csv'.")
