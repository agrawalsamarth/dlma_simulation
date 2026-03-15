#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 19:32:22 2021

@author: samarth
"""

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

df = pd.read_csv('test.csv', skiprows=19)
att_col = 'att_1'
N = len(df)
df['r_att'] = 0


for i in range(N):
    
    idx = (int)(df.iloc[i][att_col])
    
    rx = df.iloc[i]['x0'] - df.iloc[idx]['x0']
    ry = df.iloc[i]['x1'] - df.iloc[idx]['x1']
    rz = df.iloc[i]['x2'] - df.iloc[idx]['x2']
    
    df.loc[i, 'r_att'] = np.sqrt(rx**2 + ry**2 + rz**2)
    
#plt.scatter(df['x0'], df['x1'],s=1)