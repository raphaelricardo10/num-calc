import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

jacobi = np.transpose(pd.read_csv('jacobi.csv', delimiter=',', header=None))
seidel = np.transpose(pd.read_csv('seidel.csv', delimiter=',', header=None))

ij = jacobi.iloc[0]
iS = seidel.iloc[0]
x1j = jacobi.iloc[1]
x2j = jacobi.iloc[2]
x3j = jacobi.iloc[3]
x1s = seidel.iloc[1]
x2s = seidel.iloc[2]
x3s = seidel.iloc[3]

#plt.plot(ij, x1j)
#plt.plot(ij, x2j)
#plt.plot(ij, x3j)
plt.plot(iS, x1s)
plt.plot(iS, x2s)
plt.plot(iS, x3s)
plt.show()