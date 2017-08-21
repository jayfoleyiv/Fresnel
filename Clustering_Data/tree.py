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

X = np.loadtxt("/home/james/Fresnel/Clustering_Data/TEST_SET_3200K_Stop.txt", usecols=(0,1,2,3,4))
Y =  np.loadtxt("/home/james/Fresnel/Clustering_Data/TEST_SET_3200K_Stop.txt", usecols=(5,6))


loadtime = time.time()
print("load time")
print(loadtime - st)

# Fit regression model
rt = time.time()
regr = DecisionTreeRegressor()
regr.fit(X, Y)
print("time to do regression")
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


#pt = time.time()
#y_1 = regr.predict(x_1)
#print("time to predict")
#print(time.time() - pt)

pdt = time.time()
with open("regression2.pkl", 'wb') as f:
    pickle.dump(regr, f, pickle.HIGHEST_PROTOCOL)

print("pickle dump")
print(time.time()-pdt)

exptree = time.time()
export_graphviz(regr, out_file='rtree2.text')
with open("rtree2.txt", "w") as f:
    f = export_graphviz(regr, out_file=f)

print("time to export tree")
print(time.time()-exptree)

"""

#output results
print("R Squared")
print(regr.score(x_1,y_2))


ervector = np.fabs(y_1 - y_2)/y_2
vecte = ervector[0][0]
vectd = ervector[0][1]
for l in range(1, len(ervector)):
    vecte = np.vstack((vecte, ervector[l][0]))
    vectd = np.vstack((vectd, ervector[l][1]))

er =  np.linalg.norm(vecte) + np.linalg.norm(vectd)

print("Error")
print(er)

elapsed_time = time.time()-st
print("Elapsed time: %.2fs" % elapsed_time)

export_graphviz(regr, out_file='rtree2.text')
with open("rtree2.txt", "w") as f:
    f = export_graphviz(regr, out_file=f)

text_file = open("/home/james/rtree_test.txt", "a")
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
"""
print("total time")
print(time.time()-st)

