#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 14:30:30 2017

@author: jtsatsaros2018
"""

##!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 13:25:36 2017

@author: jtsatsaros2018
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import export_graphviz
import time as time
import pickle as pickle

st = time.time()

# Create a dataset

loadtime = time.time()
print("load time")
print(loadtime - st)

# Fit regression model
rt = time.time()
print(time.time()-rt)

# Predict
"""
a = list(np.arange(0.03, 0.33 ,0.03))
b = list(np.arange(6,15))
c = list(np.arange(0.19,0.28,0.01))
d = list(np.arange(0.13,0.22,0.01))
e = list(np.arange(1225,1450,25))

x_1 = np.array([0,5,0.18,0.12,1200])
for i in range(8):
    x_1 = np.vstack((x_1,[a[i],b[i],c[i],d[i],e[i]]))
"""
#x_1 = np.array([0.0,8,0.18,0.16,1200.0])


pt = time.time()

print("time to predict")
print(time.time() - pt)

x = np.array([[1],[2],[3],[4]])
y= np.array([[3],[6],[9],[12]])

regr = DecisionTreeRegressor()
regr.fit(x,y)

x_1 = np.array([[4],[5],[6],[7]])
y_1 = regr.predict(x_1)

print(y_1)

pdt = time.time()
with open("regtest.pkl", 'wb') as f:
    pickle.dump(regr, f, pickle.HIGHEST_PROTOCOL)

print("pickle dump")
print(time.time()-pt)



#output results
print("R Squared")
print("total time")
print(time.time()-st)

