import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns

sns.set()

colors = sns.color_palette().as_hex()


def my_yticks(y, pos):
    if y <= 0:
        return '$0$'
    exponent = int(np.log10(y))
    value = y / float(10 ** exponent);
    return '${{ %1.1f \mathrm{e} {%2d} }}$' % (value, exponent)


# Name of output file to save the plot to
outputFileName = 'compute_response_layers_runtime_bar_plot_greater128.png'

# Name of input files with performance data, etc.
inputFileNames = [
    '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layer.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_leftcorner.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_top.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_top_mid.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_laplacian_locality.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_laplacian_locality_unconditional.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_laplacian_locality_unconditional_opt.csv',
    # '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_Dyy_laplacian_locality_unconditional_opt_flops.csv'
    '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_switch_Dyy.csv',
    '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_unconditional.csv',
    '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_switch_Dyy_unconditional.csv',
    '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_switch_Dyy_unconditional_opt.csv',
    '../../benchmarking_files/valid_benchmarking/2020_05_23_07_36_laurin/compute_response_layers_unconditional_strided.csv'
]

# Name of labels
plotLabels = [
    'Base',
    # 'Left Corner',
    # 'Top',
    # 'Top Mid',
    # 'Dyy',
    # 'Locality',
    # 'Padded',
    # 'P Opt',
    # 'P Opt Flops'
    'Dyy',
    'padded',
    'Dyy p',
    'Dyy p opt',
    'strided'
]

# Getting current axis
ax = plt.gca()

# Initializing plot title
plt.title('Unpadded Runtime Optimizations$', x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold',
          loc='center')

# Initializing plot axis labels
# plt.xlabel('optimization versions', fontsize=10)
yl = plt.ylabel('[cycles]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

# Initializing and setting axis ticks
ticks_x = [i for i, _ in enumerate(inputFileNames)]
ticks_y = [0, 10000000, 25000000, 50000000, 75000000]

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

    print("Cycles:")
    print(avg_cycles)

    # Plotting  flops per cycles performance
    plt.bar(i, avg_cycles, label=plotLabels[i])

plt.xticks(ticks_x, plotLabels)
plt.show()
# Saving plot to file
#plt.savefig(outputFileName, dpi=300)
