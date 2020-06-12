import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns

sns.set()

colors = sns.color_palette().as_hex()


# def my_xticks(y, pos):
#     if y <= 0:
#         return '$0$'
#     exponent = int(np.log10(y))
#     value = y / float(10 ** exponent);
#     return '${{ %1.1f \mathrm{e} {%2d} }}$' % (value, exponent)
def my_xticks(x, pos):
    if x <= 0:
        return '$0$'
    value = x / int(10**6);
    return '${{ {%2d} }}$' % (value)


# Name of output file to save the plot to
outputFileName = 'response_layers_runtime_general_cache_opt.pdf'

# Name of input files with performance data, etc.
inputFileNames = [
    '../../../benchmarking_files/valid_benchmarking/2020_05_20_23_57_sebastian/compute_response_layer.csv',
    '../../../benchmarking_files/valid_benchmarking/2020_05_20_23_57_sebastian/compute_response_layers_precompute.csv',
    '../../../benchmarking_files/valid_benchmarking/2020_05_20_23_57_sebastian/compute_response_layers_blocking.csv',
    '../../../benchmarking_files/valid_benchmarking/2020_05_20_23_57_sebastian/compute_response_layers_at_once.csv',
    ]

# Name of labels
plotLabels = [
    '$\mathtt{ base }$',
    #'$\mathtt{ precomputations }$',
    '$\mathtt{ general }$',
    '$\mathtt{ blocking }$',
    '$\mathtt{ at\ once }$'
]

plotLabels.reverse()


# Getting current axis
ax = plt.gca()

# Initializing plot title
plt.title('Response Layers - General and Cache Optimizations', x=-0.18, y=1.05, ha='left', fontsize=16,
          fontweight='bold')

# Initializing plot axis labels
plt.xlabel('[mio. cycles]', fontsize=13)
#yl = plt.ylabel('Image Size 1024', fontsize=15, ha='left')
#yl.set_rotation(0)
#ax.yaxis.set_label_coords(-0.1, 1.01)

# Initializing and setting axis ticks
ticks_y = [i for i, _ in enumerate(inputFileNames)]
# ticks_y = [0, 25000000, 50000000, 75000000, 100000000] # laurin scale
ticks_x = [0, 50000000, 100000000, 150000000, 200000000] # sebastian scale

# Setting y ticks and limits (min, max)
plt.xticks(ticks_x)
plt.xlim(ticks_x[0] - 1, ticks_x[len(ticks_x) - 1] + 1)

# Setting axis ticks formatter
ax.xaxis.set_major_formatter(ticker.FuncFormatter(my_xticks))

# Setting label size
ax.tick_params(axis='both', which='major', labelsize=12)
bar_width = 0.6
# Iterating through all input files and plotting each as a line
counter = 0
inputFileNames.reverse()
for i in range(len(inputFileNames)-1, -1, -1):
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

    print("Cycles:")
    print(avg_cycles)

    # Plotting  flops per cycles performance
    #plt.bar(i, avg_cycles, label=plotLabels[i])

    p0 = plt.barh(i, avg_cycles, bar_width, color=colors[counter])
    counter += 1

#p0 = plt.barh(2, runtime_cycles[0], width, color=colors[0])
#p1 = plt.barh(1, runtime_cycles[1], width, color=colors[3])
#p1 = plt.barh(0, runtime_cycles[2], width, color=colors[2])

plt.yticks(ticks_y, plotLabels)


ratio = 0.3
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

plt.tight_layout()

plt.gcf().subplots_adjust(left=0.2)

plt.gcf().set_size_inches(7, 2.7)


#plt.show()
# Saving plot to file
plt.savefig(outputFileName, dpi=300, figsize=(500,150))
