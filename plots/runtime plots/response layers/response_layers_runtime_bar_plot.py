import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns
sns.set()

sns.set_palette("Blues")
colors = sns.color_palette().as_hex()

def my_yticks(y,pos):
    if y <= 0:
        return '$0$'
    exponent = int(np.log10(y)) 
    value = y / float(10**exponent);
    return '${{ %1.1f \mathrm{e} {%2d} }}$' % (value, exponent)


# Name of input files with performance data, etc.
inputFileNamesNoPadded = [
    '2020_05_20_23_57_response_layers/compute_response_layer.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_precompute.csv',
    '2020_05_20_23_57_response_layers/compute_response_layers_blocking.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_at_once.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_sonic_Dyy.csv', 
    #'2020_05_20_23_57_response_layers/compute_response_layers_Dyy_leftcorner.csv', 
    #'2020_05_20_23_57_response_layers/compute_response_layers_Dyy_top.csv', 
    #'2020_05_20_23_57_response_layers/compute_response_layers_Dyy_top_mid.csv', 
    #'2020_05_20_23_57_response_layers/compute_response_layers_Dyy.csv', 
    #'2020_05_20_23_57_response_layers/compute_response_layers_Dyy_laplacian.csv',
    #'2020_05_20_23_57_response_layers/compute_response_layers_Dyy_laplacian_localityloops.csv',
]

# Name of padded input files with performance data, etc.
inputFileNamesPadded = [
    '2020_05_20_23_57_response_layers/compute_response_layer.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_unconditional.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_unconditional_strided.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_sonic_Dyy_unconditional.csv', 
    '2020_05_20_23_57_response_layers/compute_response_layers_sonic_Dyy_unconditional_opt.csv' 
]

# Name of labels
plotLabelsNoPadded = [
    '$\mathtt{ base }$',
    '$\mathtt{ precompute }$', 
    '$\mathtt{ blocking }$', 
    '$\mathtt{ all\ \ at\ \ once }$', 
    '$\mathtt{ switch }\ \ D_{yy} }$' 
]

# Name of padded labels
plotLabelsPadded = [
    '$\mathtt{ base }$',
    '$\mathtt{ padded}$', 
    '$\mathtt{ padded\ \ strided }$', 
    '$\mathtt{ padded\ \ D_{yy} }$', 
    '$\mathtt{ padded\ \ D_{yy}\ \ opt }$' 
]

#inputFileNames = inputFileNamesNoPadded
#plotLabels = plotLabelsNoPadded
#outputFileName = 'response_layers_runtime_bar_plot_no_padded.png'
#title = 'Compute Response Layer Runtime Failed Optimizations'

inputFileNames = inputFileNamesPadded
plotLabels = plotLabelsPadded
outputFileName = 'response_layers_runtime_bar_plot_padded.png'
title = 'Compute Response Layer Runtime Padded Optimizations'

# Getting current axis
ax = plt.gca()

# Initializing plot title
plt.title(title,  x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.xlabel('optimization versions', fontsize=10)
yl = plt.ylabel('[cycles]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.02)

# Initializing and setting axis ticks
ticks_x = [i for i, _ in enumerate(inputFileNames)]
ticks_y = [0, 50000000, 100000000, 150000000, 200000000]

plt.xticks(ticks_x, plotLabels)

# Setting y ticks and limits (min, max)
plt.yticks(ticks_y)
plt.ylim(ticks_y[0] - 1, ticks_y[len(ticks_y) - 1] + 1)

# Setting axis ticks formatter
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
    imageName = data[5, 0]
    width = data[5, 1]
    height = data[5, 2]
    num_interest_points = data[5, 3]
    num_flops = data[5, 4]
    avg_cycles = data[5, 5]
    min_cycles = data[5, 6]
    max_cycles = data[5, 7]
    flops_per_cycles = data[5, 8]

    print(avg_cycles)

    # Plotting  flops per cycles performance
    plt.bar(i, avg_cycles, label=plotLabels[i])

 
# Saving plot to file
plt.savefig(outputFileName, dpi=300)
