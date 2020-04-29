import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
sns.set()

filenames = ['data/mvm_O0.csv', 'data/mvm_O3_fno-tree-vectorize.csv', 'data/mvm_O3_ffast-math_march-native.csv']
labels = ['$\mathtt{-O0}$', '$\mathtt{-O3}$ $\mathtt{-fno-tree-vectorize}$', '$\mathtt{-O3}$ $\mathtt{-ffast-math}$ $\mathtt{-march=naive}$']

plt.title('Matrix-Vector Multiplication Performances')
plt.xlabel('matrix size n')
plt.ylabel('[flops/cycles]')

ticks = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
plt.xticks(ticks)
plt.xlim(0, 4200)
plt.ylim(0, 2.5)

maximum = 0.0

#for i in range(0, 3):
for i in range(2, -1, -1):
    data = np.genfromtxt(filenames[i], delimiter=',')
    ns = data[:,0]
    cycles = np.median(data[:,1:], axis=1);
    performances = (2*ns*ns) / cycles
    if (max(performances) > maximum):
        maximum = max(performances)
    plt.plot(ns, performances, label=labels[i], marker='o')

plt.legend()
plt.savefig('data/mvm.png')

print("The highest performance is: " + str(maximum) + " [flops/cycle]")
