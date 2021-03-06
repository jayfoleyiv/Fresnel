import numpy as np
from sklearn.cluster import MeanShift, estimate_bandwidth
import matplotlib.pyplot as plt
from itertools import cycle
import time as time

arr = np.loadtxt("/home/james/Fresnel/Clustering_Data/TEST_SET_320000000.txt")
st = time.time()
ms = MeanShift()
ms.fit(arr)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)
elapsed_time = time.time() - st
print("Elapsed time: %.2fs" % elapsed_time)
print("number of estimated clusters : %d" % n_clusters_)


plt.figure(1)
plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(arr[my_members, 1], arr[my_members, 6], col + '.')
    plt.plot(cluster_center[1], cluster_center[6], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()
