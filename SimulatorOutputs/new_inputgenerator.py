import pandas as pd

# Load the CSV files into DataFrames
df1 = pd.read_csv('vdm_input.csv')
# df2 = pd.read_csv('steering.csv')
df2 = pd.read_csv('steering.csv')

# Specify the common column
common_column = 'TIME'  # Replace with the name of the common column

# Merge the DataFrames based on the common column
merged_df = pd.merge(df1, df2, on=common_column)

# Save the resulting DataFrame to a new CSV file
merged_df.to_csv('input_common_sample.csv', index=False)

print("Rows with common values have been merged and saved to 'merged_file.csv'.")
