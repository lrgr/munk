import sys, os, argparse, numpy as np, logging, pandas as pd
from sklearn.externals import joblib
from copy import deepcopy
import matplotlib as mp
import matplotlib.pyplot as plt

data_dir = 'collins-roguev-sc-sp/features/'

sc_data = joblib.load(data_dir + 'sc-features-no_landmarks.pkl')
sp_data = joblib.load(data_dir + 'sp-features-no_landmarks.pkl')
X_sc, y_sc = np.array(sc_data.get('X')), sc_data.get('y')
X_sp, y_sp = np.array(sp_data.get('X')), sp_data.get('y')

sc_T_sample = X_sc[np.random.choice(np.where(y_sc)[0], 1000)]
sp_T_sample = X_sp[np.random.choice(np.where(y_sp)[0], 1000)]

sc_F_sample = X_sc[np.random.choice(np.where(~y_sc)[0], 1000)]
sp_F_sample = X_sp[np.random.choice(np.where(~y_sp)[0], 1000)]

axs = [0, 1]

import seaborn as sns

# plt.figure()
sns.jointplot(sc_T_sample[:,axs[0]], sc_T_sample[:,axs[1]])
plt.title('sc - T')
sns.jointplot(sc_F_sample[:,axs[0]], sc_F_sample[:,axs[1]])
plt.title('sc - F')

# plt.figure()
sns.jointplot(sp_T_sample[:,axs[0]], sp_T_sample[:,axs[1]])
plt.title('sp - T')
sns.jointplot(sp_F_sample[:,axs[0]], sp_F_sample[:,axs[1]])
plt.title('sp - F')

plt.show()
