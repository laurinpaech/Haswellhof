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
    return r'${{ {:2d} }}$'.format(int(y))

# Name of output file to save the plot to
outputFileName = 'compute_response_layers_runtime_plot.png'

# Name of input files with performance data, etc.
inputFileNames = [
    '../../benchmark/2020_05_20_23_57_first_try/compute_response_layer.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_precompute.csv',
    '../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_blocking.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_at_once.csv', 
    '../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_sonic_Dyy.csv', 
    #'../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_Dyy_leftcorner.csv', 
    #'../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_Dyy_top.csv', 
    #'../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_Dyy_top_mid.csv', 
    #'../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_Dyy.csv', 
    #'../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_Dyy_laplacian.csv',
    #'../../benchmark/2020_05_20_23_57_first_try/compute_response_layers_Dyy_laplacian_localityloops.csv'
]

# Name of labels
plotLabels = [
    '$\mathtt{ base }$', 
    '$\mathtt{ precompute }$',
    '$\mathtt{ blocking }$',
    '$\mathtt{ at\ once }$',
    '$\mathtt{ sonic\ Dyy }$', 
    #'$\mathtt{ Dyy\ leftcorner }$', 
    #'$\mathtt{ Dyy\ top }$', 
    #'$\mathtt{ Dyy\ top\ mid }$', 
    #'$\mathtt{ Dyy }$', 
    #'$\mathtt{ Dyy\ laplacian }$',
    #'$\mathtt{ Dyy\ laplacian locality }$'
]

# Getting current axis
ax = plt.gca()

# Initializing plot title
plt.title('Runtime of $\mathtt{ compute\_response\_layers }$ optimizations',  x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.xlabel('image size', fontsize=10)
yl = plt.ylabel('[cycles]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

# Setting x-axis to be log axis
plt.xscale('log')
#plt.yscale('log')


# Initializing and setting axis ticks
ticks_x = []
for i in range(8, 11):
    ticks_x.append(2**(i))
ticks_y = [0, 50000000, 100000000, 150000000, 200000000]

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
    imageName = data[3:6, 0]
    width = data[3:6, 1]
    height = data[3:6, 2]
    num_interest_points = data[3:6, 3]
    num_flops = data[3:6, 4]
    avg_cycles = data[3:6, 5]
    min_cycles = data[3:6, 6]
    max_cycles = data[3:6, 7]
    flops_per_cycles = data[3:6, 8]

    print("Cycles:")
    print(avg_cycles)

    # Plotting  flops per cycles performance
    plt.plot(width, avg_cycles, label=plotLabels[i], marker='o')

# Adding legend to plot
plt.legend(loc=1, prop={'size': 6})

# Saving plot to file
plt.savefig(outputFileName, dpi=300)
