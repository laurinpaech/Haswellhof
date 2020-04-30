import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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

def myticks(x,pos):
    return r'${{ {:2d} }}$'.format(x)

# Name of output file to save the plot to
outputFileName = 'performance_plot.png'

# Name of input files with performance data, etc.
inputFileNames = ['2020_04_29_16_11/create_integral_img.csv', '2020_04_29_16_11/interpolate_step.csv']

# Name of labels
plotLabels = ['$\mathtt{ create\_integral\_img }$', '$\mathtt{ interpolate\_step }$']

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
for i in range(5, 13):
    ticks_x.append(2**(i))
ticks_y = [0, 2, 4, 6]

plt.xticks(ticks_x, va='center')
plt.yticks(ticks_y)

# Setting x and y limits (min, max)
plt.xlim(ticks_x[0] / 2.0, ticks_x[len(ticks_x) - 1] * 2.0)
plt.ylim(ticks_y[0] - 1, ticks_y[len(ticks_y) - 1] + 1)

# Setting axis ticks formatter
ax.xaxis.set_major_formatter(ticker.FuncFormatter(myticks))
#ax.yaxis.set_major_formatter(ticker.FuncFormatter(myticks))

# Setting label size
ax.tick_params(axis='both', which='major', labelsize=8)

# Iterating through all input files and plotting each as a line
for i in range(0, len(inputFileNames)):

    substrings = inputFileNames[i].split('.')

    # Getting function name from file name
    functionName = substrings[0]

    # Reading csv data from input file 
    data = np.genfromtxt(inputFileNames[i], delimiter=',')

    # Getting width, height, number of interest points, average/min/max cycles, flops per cycles
    imageName = data[:, 0]
    width = data[:, 1]
    height = data[:, 2]
    num_interest_points = data[:, 3]
    num_flops = data[:, 4]
    avg_cycles = data[:, 5]
    min_cycles = data[:, 6]
    max_cycles = data[:, 7]
    flops_per_cycles = data[:, 8]

    # Plotting  flops per cycles performance
    plt.plot(width, flops_per_cycles, label=plotLabels[i], color=colors[i], marker='o')

# Adding legend to plot
plt.legend()

# Saving plot to file
plt.savefig(outputFileName, dpi=300)
