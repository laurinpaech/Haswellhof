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

def my_xticks(x,pos):
    exponent = np.log2(x)
    if exponent < 0.0:
        value = int(2**(-exponent))
        return r'$1 / {den:2d}$'.format(den=value)
    else:
        value = int(2**(exponent))
        return r'${{ {:2d} }}$'.format(value)

def my_yticks(y,pos):
    return r'${{ {:2d} }}$'.format(y)

# Name of output file to save the plot to
outputFileName = 'get_msurf_descriptor_runtime_plot.png'

# Name of input files with performance data, etc.
inputFileNames = [
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_improved.csv',
    #'../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_improved_flip.csv', 
    #'../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_improved_flip_flip.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_inlined.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_inlinedHaarWavelets.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_inlinedHaarWavelets_precheck_boundaries.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_gauss_compute_once_case.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_gauss_pecompute_haar.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_gauss_pecompute_haar_rounding.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_arrays.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/get_msurf_descriptors_haar_unroll_2_24_True_winner.csv'
]

# Name of labels
plotLabels = [
    '$\mathtt{ base }$', 
    '$\mathtt{ improved }$',
    #'$\mathtt{ improved flip }$',
    #'$\mathtt{ improved flip flip }$',
    '$\mathtt{ inlined }$',
    '$\mathtt{ joined\ haar\ wavelets }$',
    '$\mathtt{ precheck\ boundaries }$', 
    '$\mathtt{ gauss\ (cases) }$', 
    '$\mathtt{ precompute\ haar }$', 
    '$\mathtt{ rounding }$', 
    '$\mathtt{ gauss\ (arrays) }$', 
    '$\mathtt{ blocked\ 2x24\ autotune }$'
]

# Getting current axis
ax = plt.gca()

# Initializing plot title
plt.title('Runtime of $\mathtt{ get\_msurf\_descriptors }$ optimizations',  x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.xlabel('image size', fontsize=10)
yl = plt.ylabel('[cycles]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

# Setting x-axis to be log axis
plt.xscale('log')
plt.yscale('log')


# Initializing and setting axis ticks
ticks_x = []
for i in range(6, 13):
    ticks_x.append(2**(i))
ticks_y = [10000, 100000, 1000000]

plt.xticks(ticks_x, va='center')
plt.yticks(ticks_y)

# Setting x and y limits (min, max)
plt.xlim(ticks_x[0] / 2.0, ticks_x[len(ticks_x) - 1] * 2.0)
plt.ylim(ticks_y[0] - 1, ticks_y[len(ticks_y) - 1] + 1)
#plt.ylim(ticks_y[0] / 2.0, ticks_y[len(ticks_y) - 1] * 2.0)

# Setting axis ticks formatter
ax.xaxis.set_major_formatter(ticker.FuncFormatter(my_xticks))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(my_yticks))

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
    imageName = data[1:, 0]
    width = data[1:, 1]
    height = data[1:, 2]
    num_interest_points = data[1:, 3]
    num_flops = data[1:, 4]
    avg_cycles = data[1:, 5]
    min_cycles = data[1:, 6]
    max_cycles = data[1:, 7]
    flops_per_cycles = data[1:, 8]

    print("Cycles:")
    print(avg_cycles)

    # Plotting  flops per cycles performance
    plt.plot(width, avg_cycles, label=plotLabels[i], marker='o')

# Adding legend to plot
plt.legend(loc=1, prop={'size': 6})

# Saving plot to file
plt.savefig(outputFileName, dpi=300)
