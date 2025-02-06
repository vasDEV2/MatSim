from matplotlib import pyplot as plt
import pandas as pd
import math

# df1 = pd.read_csv('new.csv')
# df2 = pd.read_csv('/home/vasudevan/VDM/trainigbkc/training_data6.csv')
# df3 = pd.read_csv('/home/vasudevan/VDM/trainigbkc/training_data3.csv')
# merg_df = pd.concat([df1,df2])
# merg_df.to_csv('/home/vasudevan/VDM/trainigbkc/merg_df.csv')
df1 = pd.read_csv('input_common_sample.csv')
df2 = pd.read_csv('predict_compare.csv')
# print(df1['FXFL'])
# print(df1['lamdaFL'])
# plt.scatter(df1['lamdaFL'].to_numpy(),df1['fxFL'].to_numpy())

# plt.scatter(-df1['tauFR'].to_numpy(),df1['FYFR'].to_numpy())

# indices = df1[df1['lamdaFL'] < -0.02].index
# # 
# print(indices)
# # 
# df1.loc[indices, 'FXFL'] = df1.loc[indices, 'FXFL'].apply(
#     lambda x: None if x > -1500 else x
# )
# # 
# # 
# indices2 = df1[df1['lamdaFL'] > 0.02].index
# # 
# print(indices)
# # 
# df1.loc[indices2, 'FXFL'] = df1.loc[indices2, 'FXFL'].apply(
#     lambda x: None if x > -1500 else x
# )
# 
# df1.dropna(subset=['FYFR'], inplace=True)

# df1['tauFR'] = -df1['tauFR']

# df1.to_csv('train_tau_updated.csv')

fig, axs = plt.subplots(2,1)
# axs[0].scatter(df2['Ynew'].to_numpy(),df2['Xnew'].to_numpy(),label = 'prediction')
# axs[0].scatter(df2['Y'].to_numpy(),df2['X'].to_numpy(),label = 'actual')
# plt.scatter(df2['SLIPFL'].to_numpy(),df1['FXFL'].to_numpy())
# fig, axs = plt.subplots(2,1)
# axs.scatter(df2['Y'].to_numpy(),df2['X'].to_numpy())
axs[0].plot((df2['AX_NEW']*0.5).to_numpy(), label = 'Predction')
axs[0].plot((df2['AX']).to_numpy(), label = 'actual')
axs[1].plot((df2['AY_NEW']).to_numpy(), label = 'prediction')
axs[1].plot((df2['AY']).to_numpy(), label = 'actual')
# axs[0].plot(df2['AX'])
axs[0].legend()
# axs.axis('off')

# Set the figure background to transparent
# fig.patch.set_alpha(0)  # The overall background (outside the axes)
# axs.patch.set_alpha(0)
# az = (df1['fxFL']/df2['FFL']).mean()
# print(az)

plt.show()
