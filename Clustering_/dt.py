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

# Create a dataset
arr = np.loadtxt("/Users/jtsatsaros2018/Documents/Fresnel/Fresnel/Clustering_Data/TEST_SET_100000.txt")
X = np.array([arr[0][:5]])
Y = np.array([arr[0][5:]])
for n in range(1, len(arr)):
    X = np.vstack((X, arr[n][:5]))
    Y = np.vstack((Y, arr[n][5:]))

newarr = np.loadtxt("/Users/jtsatsaros2018/Documents/Fresnel/TEST_SET_320000000.txt")
x_1 = np.array([newarr[4][:5]])
y_2 = np.array([newarr[4][5:]])
for n in range(5, len(newarr)):
    if (newarr[n][:5] not in X ):
        x_1 = np.vstack((x_1, newarr[n][:5]))
        y_2 = np.vstack((y_2, newarr[n][5:]))


# Fit regression model
regr = DecisionTreeRegressor()
regr.fit(X, Y)


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
y_1 = regr.predict(x_1)

print(regr.score(x_1,y_2))
export_graphviz(regr, out_file='tree.dot')



with open("tree.txt", "w") as f:
    f = export_graphviz(regr, out_file=f)

"""
#output the results
text_file = open("/Users/jtsatsaros2018/Documents/test2", "a")
for k in range(len(y_1)):
    for l in range(2):
        n = y_1[k][l]
        text_file.write(str(n))
        text_file.write("\t")
    for j in range(2):
        n = y_2[k][j]
        text_file.write(str(n))
        text_file.write("\t")
    text_file.write("\n")    
text_file.close()
"""

# Plot the results
plt.figure()
s = 50
s = 25
plt.scatter(y_1[:,0],y_1[:,1], c="navy", s=s,
            edgecolor="black", label="prediction")
plt.scatter(y_2[:,0],y_2[:,1], c="red", s=s,
            edgecolor="black", label="dataset")

#plt.xlim([-6, 6])
#plt.ylim([-6, 6])
plt.xlabel("target 1")
plt.ylabel("target 2")
plt.title("Multi-output Decision Tree Regression")
plt.legend(loc="best")
plt.show()
