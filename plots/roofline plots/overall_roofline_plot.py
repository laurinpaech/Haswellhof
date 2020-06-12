import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from matplotlib.patches import Patch


import seaborn as sns
sns.set()
colors = sns.color_palette().as_hex()

def add_boundaries(ax, plot_min, plot_max, color):
    scalar_bound_x = np.array([plot_min, plot_max])
    scalar_bound_y = np.array([4, 4])

    avx2_bound_x = np.array([plot_min, plot_max])
    avx2_bound_y = np.array([32, 32])

    mem_bound_x = np.array([plot_min, plot_max])
    mem_bound_y = 5.62 * mem_bound_x

    cache_L3_bound_x = np.array([plot_min, plot_max])
    cache_L3_bound_y = (19.32/2.7) * cache_L3_bound_x

    cache_L2_bound_x = np.array([plot_min, plot_max])
    cache_L2_bound_y = (23.48/2.7) * cache_L2_bound_x

    cache_L1_bound_x = np.array([plot_min, plot_max])
    cache_L1_bound_y = (57.8/2.7) * cache_L1_bound_x


    ax.plot(scalar_bound_x, scalar_bound_y, c=color, zorder=1, linewidth=0.5)
    ax.plot(avx2_bound_x, avx2_bound_y, c=color, zorder=1, linewidth=0.5)
    ax.plot(mem_bound_x, mem_bound_y, c=color, zorder=1, linewidth=0.5)
    ax.plot(cache_L3_bound_x, cache_L3_bound_y, c=color, zorder=1, linewidth=0.5)    
    ax.plot(cache_L2_bound_x, cache_L2_bound_y, c=color, zorder=1, linewidth=0.5)    
    ax.plot(cache_L1_bound_x, cache_L1_bound_y, c=color, zorder=1, linewidth=0.5)    

    scalar_label_pos = (1.5 * 1.1, 4.0 * 1.01)
    avx2_label_pos = (1.5 * 1.3, 32.0 * 1.01)
    mem_label_pos = (0.95, 5.62 * 1.01)
    cache_L3_label_pos = (0.95, 7.16 * 1.01)
    cache_L2_label_pos = (0.95, 8.70 * 1.01)
    cache_L1_label_pos = (0.95, 21.41 * 1.01)


    ax.annotate('$\\pi_{scalar} = 4$', xy=scalar_label_pos, va='bottom', fontsize=7)
    ax.annotate('$\\pi_{AVX2} = 32$', xy=avx2_label_pos, va='bottom', fontsize=7)
    #ax.annotate('$\\mathtt{P} \\leq \\beta \\mathtt{I} = 32 \\mathtt{I}$', xy=mem_label_pos, ha='center', va='center', fontsize=8, rotation=45)
    ax.annotate('$\\mathtt{DRAM}:\ 5.62$', xy=mem_label_pos, ha='center', va='center', rotation=34, fontsize=7)
    ax.annotate('$\\mathtt{L3}:\ 7.16$', xy=cache_L3_label_pos, ha='center', va='center', rotation=34, fontsize=7)
    ax.annotate('$\\mathtt{L2}:\ 8.70$', xy=cache_L2_label_pos, ha='center', va='center', rotation=34, fontsize=7)
    ax.annotate('$\\mathtt{L1}:\ 21.41$', xy=cache_L1_label_pos, ha='center', va='center', rotation=34, fontsize=7)



def add_point(ax, x, y, label, label_offset, c, m):
    handle = ax.scatter(x, y, color=c, zorder=2, marker=m)
    label_pos = (x * (1.0 + label_offset[0]), y * (1.0 + label_offset[1]))
    ax.annotate(label, xy=label_pos, va='center', fontsize=8)
    return handle

def myticks(x,pos):
    if x == 0: return "$0$"
    return '$ %1.1f $' % (x)


'''
def myticks(x,pos):
    if x == 0: return "$0$"
    exponent = np.log2(x)
    #coeff = x/2**exponent
    if exponent < 0.0:
        value = int(2**(-exponent))
        s = '$\\frac{ 1 }{ %d }$' % (value)
        return s
    else:
        value = int(2**(exponent))
        s = '$ %d $' % (value)
        return s
'''

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
'''
def myticks(x,pos):
    if x == 0: return "$0$"
    exponent = int(np.log10(x))
    coeff = x/(10**exponent)
    if (coeff == 1.0) or (coeff == 0.5):
        return r"${{ {:2.1f} }}$".format(x)
    else:
        return ""
'''


ax = plt.gca()

plt.title('Roofline Plot - Descriptors and Response Layers', x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')
plt.xlabel('Operational Intensity $\\mathtt{I}$ [flops/byte]', fontsize=11)
yl = plt.ylabel('Performance $\\mathtt{P}$ [flops/cycle]', fontsize=11, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

plt.xscale('log')
plt.yscale('log')

#ticks_x = [2**(-4), 2**(-3),2**(-2), 2**(-1), 2**(0), 2**(1), 2**(2), 2**(3), 2**(4)]
#ticks_y = [2**(0), 2**(1), 2**(2), 2**(3), 2**(4), 2**(5)]
#ticks_x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0]
#ticks_y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0]

#plt.xticks(ticks_x, va='center')
#plt.yticks(ticks_y)

plt.xlim(10**(-2), 5)
plt.ylim(2**(-4), 64)

ax.xaxis.set_major_formatter(ticker.FuncFormatter(myticks))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(myticks))

#ax.tick_params(axis='both', which='major', labelsize=8)

# Adding boundaries from 3 (a)
add_boundaries(ax, 10**(-6), 10**(6), sns.xkcd_rgb['grey'])

# Adding intersecition points of boundaries from 3 (a)
#add_point(ax, 5.0/32.0, 5.0, '', (0., 0.0), 'grey', 'o')
#add_point(ax, 20.0/32.0, 20.0, '', (0., 0.0), 'grey', 'o')


data = pd.read_csv('data_fixed.csv')

print(data)

runtime = data['runtime']

# Multiplying times number of interest points
for i in range(0, data.shape[0]):
    if (i < 9):
        runtime[i] *= 2846


work = data['flops']
data_movement = data['memory']

labels = [
    # "$\mathtt{ (b)\ iimage }$",
    # "$\mathtt{ (b)\ response }$",
    # "$\mathtt{ (b)\ interest }$",
    # "$\mathtt{ (b)\ descriptors }$",
    # "$\mathtt{ (o)\ iimage }$",
    # "$\mathtt{ (o)\ response }$",
    # "$\mathtt{ (o)\ interest }$",
    # "$\mathtt{ (o)\ descriptors }$",
    # "$\mathtt{ (op)\ iimage }$",
    # "$\mathtt{ (op)\ response }$",
    # "$\mathtt{ (op)\ interest }$",
    # "$\mathtt{ (op)\ descriptors }$"
    '$\mathtt{(d)\ baseline}$',
    '$\mathtt{(d)\ inlined}$',
    '$\mathtt{(d)\ haarXY}$',#2
    '$\mathtt{(d)\ reduced\ exp()}$', #3
    '$\mathtt{(d)\ precompute\ haar}$', #4
    '$\mathtt{(d)\ rounding}$',#5
    '$\mathtt{(d)\ unroll\ 2x24}$',#6
    '$\mathtt{(d)\ avx2}$',#7
    '$\mathtt{(d)\ avx2\ 2x24}$',
    '$\mathtt{(rl)\ baseline}$',
    '$\mathtt{(rl)\ D_{yy}}$',
    '$\mathtt{(rl)\ padding\ D_{yy}}$'
]

label_offset = [
    # (0.1, 0.0),
    # (0.1, 0.0),
    # (0.1, 0.0),
    # (0.1, 0.0),
    # (0.1, 0.0),
    # (0.1, 0.0),#5
    # (-0.45, 0.12),#6
    # (0.1, -0.02),#7
    # (0.1, 0.0),
    # (-0.57, 0),
    # (-0.36, 0.15),
    # (-0.62, -0.05)

    (0.1, -0.05),
    (0.1, 0.0),
    (0.1, -0.03),#2
    (0.1, 0.0), #3
    (0.1, -0.02), #4
    (0.1, -0.06),  # 5
    (-0.67, 0.12),  # 6
    (0.1, -0.04),  # 7
    (0.1, 0.02),
    (-0.67, -0.04),
    (-0.43, 0.15),
    (-0.68, -0.11)
]


markers = ['o', '^', '*']
color = 0
marker = 0

handles = []
legend_labels = ['$\mathtt{(d)\ \ descriptors}$', '$\mathtt{(d)\ \ descriptors\ \ avx2}$', '$\mathtt{(rl)\ \ response\ \ layer}$']
for i in range(0, data.shape[0]):
    intensity = work[i] / data_movement[i]
    performance = work[i] / runtime[i]
    #add_point(ax, intensity, performance, labels[i], label_offset[i], colors[i % 4], markers[i // 4])
    if(i == 7 or i == 8):
        color = 2
    elif(i > 8):
        color = 1
        marker = 1
    h = add_point(ax, intensity, performance, labels[i], label_offset[i], colors[color], markers[marker])
    if(i == 0 or i == 7 or i == 9):
        handles.append(h)
    print('(%f, %f)' % (intensity, performance))

# legend_elements = [Patch(colors[0], marker=markers[0], lw=4, label='descriptors'),
#                    Patch([0],[1], marker=markers[0], color=colors[2], label='descriptors SIMD',
#                           markerfacecolor=colors[2], markersize=15),
#                    Patch([0],[2], marker=markers[1], color=colors[1], label='compute response layer',
#                           markerfacecolor=colors[1], markersize=15)]
#
#
plt.legend(handles= handles, labels=legend_labels, loc='lower right', fontsize=10)

plt.savefig('overall_roofline_plot.pdf')

