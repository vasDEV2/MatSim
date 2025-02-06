import pandas as pd

df = pd.read_csv("vdm_input.csv")

df['Fz_1'] = 0.022129
df['Fz_2'] = 0.022129
df['Fz_3'] = 0.022129
df['Fz_4'] = 0.022129

df.to_csv('yo.csv')