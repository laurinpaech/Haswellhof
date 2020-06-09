import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns
sns.set()

colors = sns.color_palette().as_hex()

def my_yticks(y,pos):
    if y <= 0:
        return '$0$'
    exponent = int(np.log10(y)) 
    value = y / float(10**exponent);
    return '${{ %1.1f \mathrm{e} {%2d} }}$' % (value, exponent)

# Getting current axis
ax = plt.gca()

outputFileName = 'overall_stacked_bar_plot.png'

# Initializing plot title
plt.title('SURF Program Runtime Comparison',  x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.xlabel('', fontsize=10)
yl = plt.ylabel('[cycles]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

# Setting axis ticks formatter
ax.yaxis.set_major_formatter(ticker.FuncFormatter(my_yticks))

N = 3

# Name of labels
plotLabels = [
    '$\mathtt{ compute\_integral\_image }$',
    '$\mathtt{ compute\_response\_layers }$',
    '$\mathtt{ get\_interest\_points }$',
    '$\mathtt{ get\_msurf\_descriptors }$'
]

'''
base:
4257082, 4244716, 9363333
203633614, 203087292, 212466254
26371938, 26269584, 29461552
610610269, 609645028, 644719928
optimized:
2764978, 2735236, 7650836
212745057, 212102712, 224582853
26383849, 26255828, 28194101
30229192, 29969980, 32834080
optimized padded:
4591420, 4524704, 11597536
44896938, 44594472, 52467567
26421902, 26274404, 32644300
28272216, 27939780, 37983852
optimized (v2):
3524510, 3501016, 9100796
211219059, 210404116, 224450048
26390002, 26228000, 28061004
30725934, 30501196, 33853412
'''

# TODO: Redo benchmark for more accurate results
# base, optimized (unpadded), optimized (padded)
#compute_integral_image_cycles = [4096364, 2772316, 3933956]
#compute_response_layers_cycles = [135464910, 132424959, 44925334]
#get_interest_points_cycles = [30326909, 30393429, 26568530]
#get_msurf_descriptors_cycles = [569110706, 44537542, 53325326]

compute_integral_image_cycles = [4257082, 2764978, 4591420]
compute_response_layers_cycles = [203633614, 212745057, 44896938]
get_interest_points_cycles = [26371938, 26383849, 26421902]
get_msurf_descriptors_cycles = [610610269, 30229192, 28272216]

ind = np.arange(N)    # the x locations for the groups
width = 0.4           # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, get_msurf_descriptors_cycles, width, color=colors[3])
p2 = plt.bar(ind, get_interest_points_cycles, width, bottom=get_msurf_descriptors_cycles, color=colors[2])
p3 = plt.bar(ind, compute_response_layers_cycles, width, bottom=np.array(get_msurf_descriptors_cycles)+np.array(get_interest_points_cycles), color=colors[1])
p4 = plt.bar(ind, compute_integral_image_cycles, width, bottom=np.array(get_msurf_descriptors_cycles)+np.array(get_interest_points_cycles)+np.array(compute_response_layers_cycles), color=colors[0])

plt.xticks(ind, ('$\mathtt{ base }$', '$\mathtt{ optimized }$\n$\mathtt{ (not\ \ padded) }$', '$\mathtt{ optimized }$\n$\mathtt{ (padded) }$'))
plt.yticks(np.arange(0, 800000001, 100000000))
plt.legend((p4[0], p3[0], p2[0], p1[0]), plotLabels)
#plt.legend((p1[0], p2[0], p3[0]), ('1', '2', '3'))

# Saving plot to file
plt.savefig(outputFileName, dpi=300)
