import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import figaspect

import seaborn as sns
sns.set()

colors = sns.color_palette().as_hex()

'''
def my_xticks(x,pos):
    if x <= 0:
        return '$0$'
    exponent = int(np.log10(x)) 
    value = x / float(10**exponent);
    return '${{ %1.1f \mathrm{e} {%2d} }}$' % (value, exponent)
'''
def my_xticks(x,pos):
    if x <= 0:
        return '$0$'
    value = x / int(10**6);
    return '${{ {%2d} }}$' % (value)


# Getting current axis
ax = plt.gca()

outputFileName = 'implementation_runtime_bar_plot.pdf'

# Initializing plot title
plt.title('Runtime - Different SURF Implementations',  x=-0.22, y=1.05, ha='left', fontsize=16, fontweight='bold')

# Initializing plot axis labels
plt.ylabel('', fontsize=10, labelpad=100)
xl = plt.xlabel('[mio. cycles]', fontsize=10)
#yl.set_rotation(0)
#ax.xaxis.set_label_coords(-0.1, 1.00)

#plt.xscale('log')
#plt.xticks([5*10**7, 10**8, 5*10**8])
#plt.xlim(10**7, 10**9)

# Setting axis ticks formatter
ax.xaxis.set_major_formatter(ticker.FuncFormatter(my_xticks))

N = 2

# Name of labels
plotLabels = [
    'OpenSURF [2]',
    'Herbert Bay [1]',
    'Ours'
]

runtime_cycles = [764998891, 312253445, 108804305]


ind = np.arange(N, -1, -1)    # the x locations for the groups
print(ind)
width = 0.6           # the width of the bars: can also be len(x) sequence

p0 = plt.barh(2, runtime_cycles[0], width, color=colors[0])
p1 = plt.barh(1, runtime_cycles[1], width, color=colors[3])
p1 = plt.barh(0, runtime_cycles[2], width, color=colors[2])

plt.yticks(ind, ('OpenSURF [2]', 'Herbert Bay [1]', 'Ours'))

ratio = 0.5
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

plt.tight_layout()

plt.gcf().subplots_adjust(left=0.2)

plt.gcf().set_size_inches(7, 3.6)

# Saving plot to file
plt.savefig(outputFileName, dpi=300, figsize=(500, 200))
