import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
sns.set()

colors = sns.color_palette().as_hex()

'''
def myticks(x,pos):
    if x == 0: return "$0$"
    exponent = np.log2(x)
    #coeff = x/2**exponent
    if exponent < 0.0:
        value = int(2**(-exponent))
        return r'$1 / {den:2d}$'.format(den=value)
    else:
        value = int(2**(exponent))
        return r'${{ {:2d} }}$'.format(value)
'''

# Name of output file to save the plot to
outputFileName = 'performance_plot.png'

# Name of input files with performance data, etc.
inputFileNames = ['img1-get_integral_image', 'img2-get_interest_points']

# Name of labels
plotLabels = ['$\mathtt{ get_integral_image }$, img1', '$\mathtt{ get_interest_points }$, img2']

# Getting current axis
ax = plt.gca()

# Initializing plot title
plt.title('Performance Plot',  x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.xlabel('image size', fontsize=10)
yl = plt.ylabel('[flops/cycles]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

# Setting x-axis to be log axis
plt.xscale('log')

# Initializing and setting axis ticks
ticks_x = []
for i in range(4, 13):
    ticks_x.append(2**(i))
ticks_y = [0, 2, 4, 6, 8, 10]

plt.xticks(ticks_x, va='center')
plt.yticks(ticks_y)

# Setting x and y limits (min, max)
plt.xlim(ticks_x[0] / 2.0, ticks_x[len(ticks_x) - 1] * 2.0)
plt.ylim(ticks_y[0] - 1, ticks_y[len(ticks_y) - 1] + 1)

# Setting axis ticks formatter
#ax.xaxis.set_major_formatter(ticker.FuncFormatter(myticks))
#ax.yaxis.set_major_formatter(ticker.FuncFormatter(myticks))

# Setting label size
ax.tick_params(axis='both', which='major', labelsize=8)

# Iterating through all input files and plotting each as a line
for i in range(0, len(inputFileNames)):

    substrings = inputFileNames[i].split('-')

    # Getting image name and function name from file name
    imageName = substrings[0]
    functionName = substrings[1]

    # Reading csv data from input file 
    data = np.genfromtxt(inputFileNames[i], delimiter=',')

    # Getting width, height, number of interest points, average/min/max cycles, flops per cycles
    width = data[:, 0]
    height = data[:, 1]
    num_interest_points = data[:, 2]
    avg_cycles = data[:, 3]
    min_cycles = data[:, 4]
    max_cycles = data[:, 5]
    flops_per_cycles = data[:, 6]
    flops = data[:, 7]

    # Plotting  flops per cycles performance
    plt.plot(width, flops_per_cycle, label=plotLabels[i], color=colors[i], marker='o')

# Adding legend to plot
plt.legend()

# Saving plot to file
plt.savefig(outputFileName, dpi=300)
