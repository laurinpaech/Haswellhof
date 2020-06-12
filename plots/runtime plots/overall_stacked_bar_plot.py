import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import figaspect

import seaborn as sns
sns.set()

colors = sns.color_palette().as_hex()

def my_xticks(x,pos):
    if x <= 0:
        return '$0$'
    #exponent = int(np.log10(y)) 
    value = x / int(10**6);
    return '${{ {%2d} }}$' % (value)

# Getting current axis
ax = plt.gca()

outputFileName = 'overall_stacked_bar_plot.pdf'

# Initializing plot title
plt.title('SURF Runtime Comparisons',  x=-0.22, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.ylabel('', fontsize=10, labelpad=100)
xl = plt.xlabel('[mio. cycles]', fontsize=10)
#yl.set_rotation(0)
#ax.xaxis.set_label_coords(-0.1, 1.00)

# Setting axis ticks formatter
ax.xaxis.set_major_formatter(ticker.FuncFormatter(my_xticks))

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

# compute_integral_image_cycles = [4257082, 2764978, 4591420]
# compute_response_layers_cycles = [203633614, 212745057, 44896938]
# get_interest_points_cycles = [26371938, 26383849, 26421902]
# get_msurf_descriptors_cycles = [610610269, 30229192, 28272216]

compute_integral_image_cycles = [4058780, 2785559, 4591420]
compute_response_layers_cycles = [135342608, 129156728, 45818320]
get_interest_points_cycles = [27068836, 27068836, 27068836]
get_msurf_descriptors_cycles = [568098598, 29552864, 27964796]
#overall_cycles = [743217661, 194846528, 108804305]

opensurf_cycles = [764998891]
herbertbay_cycles = [312253445]


ind = np.arange(N, 0, -1)    # the x locations for the groups
width = 0.6           # the width of the bars: can also be len(x) sequence

#p0 = plt.barh(ind, overall_cycles, width, color='grey')
p1 = plt.barh(ind, get_msurf_descriptors_cycles, width, color=colors[3])
p2 = plt.barh(ind, get_interest_points_cycles, width, left=get_msurf_descriptors_cycles, color=colors[2])
p3 = plt.barh(ind, compute_response_layers_cycles, width, left=np.array(get_msurf_descriptors_cycles)+np.array(get_interest_points_cycles), color=colors[1])
p4 = plt.barh(ind, compute_integral_image_cycles, width, left=np.array(get_msurf_descriptors_cycles)+np.array(get_interest_points_cycles)+np.array(compute_response_layers_cycles), color=colors[0])

plt.yticks(ind, ('$\mathtt{ baseline }$', '$\mathtt{ optimized }$\n$\mathtt{ (not\ \ padded) }$', '$\mathtt{ optimized }$\n$\mathtt{ (padded) }$'))
plt.xticks(np.arange(0, 800000001, 100000000))
plt.legend((p4[0], p3[0], p2[0], p1[0]), plotLabels)

ratio = 0.5
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

plt.tight_layout()

plt.gcf().subplots_adjust(left=0.2)

plt.gcf().set_size_inches(7, 3.6)

# Saving plot to file
plt.savefig(outputFileName, dpi=300, figsize=(500, 200))
