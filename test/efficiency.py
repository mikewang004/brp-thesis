import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("data-133-144-eff.txt", skiprows = 148, usecols=[1,2,3,4,5])
#print(data)

plt.scatter(data[:, 1], data[:, 2])
#plt.xlim(0.4, 1.2)
#plt.ylim(0.5, 1.5)
#print(data[:, 2])