import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns
sns.set()
colors = sns.color_palette().as_hex()

def add_boundaries(ax, plot_min, plot_max, color):
    scalar_bound_x = np.array([plot_min, plot_max])
    scalar_bound_y = np.array([5.0, 5.0])

    simd_bound_x = np.array([plot_min, plot_max])
    simd_bound_y = np.array([20.0, 20.0])

    mem_bound_x = np.array([plot_min, plot_max])
    mem_bound_y = 32.0 * mem_bound_x

    ax.plot(scalar_bound_x, scalar_bound_y, c=color, zorder=1)
    ax.plot(simd_bound_x, simd_bound_y, c=color, zorder=1)
    ax.plot(mem_bound_x, mem_bound_y, c=color, zorder=1)

    scalar_label_pos = (8.0 * 1.1, 5.0 * 1.1)
    simd_label_pos = (8.0 * 1.1, 20.0 * 1.1)
    mem_label_pos = (0.9, 32.0 * 1.1)

    ax.annotate('$\\mathtt{P} \\leq \\pi_{scalar} = 5$', xy=scalar_label_pos, va='bottom', fontsize=8)
    ax.annotate('$\\mathtt{P} \\leq \\pi_{SIMD} = 20$', xy=simd_label_pos, va='bottom', fontsize=8)
    ax.annotate('$\\mathtt{P} \\leq \\beta \\mathtt{I} = 32 \\mathtt{I}$', xy=mem_label_pos, ha='center', va='center', fontsize=8, rotation=45)



def add_point(ax, x, y, label, label_offset, color):
    handle = ax.scatter(x, y, c=color, zorder=2)
    label_pos = (x * (1.0 + label_offset[0]), y * (1.0 + label_offset[1]))
    ax.annotate(label, xy=label_pos, va='center', fontsize=8)
    return handle

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

plt.title('Roofline Plot for $\\mathtt{computation}$ Functions', x=-0.1, y=1.05, ha='left', fontsize=16, fontweight='bold')
plt.xlabel('Operational Intensity $\\mathtt{I}$ [flops/byte]', fontsize=10)
yl = plt.ylabel('Performance $\\mathtt{P}$ [flops/cycle]', fontsize=10, ha='left')
yl.set_rotation(0)
ax.yaxis.set_label_coords(-0.1, 1.01)

plt.xscale('log')
plt.yscale('log')

ticks_x = [2**(-4), 2**(-3),2**(-2), 2**(-1), 2**(0), 2**(1), 2**(2), 2**(3), 2**(4)]
ticks_y = [2**(0), 2**(1), 2**(2), 2**(3), 2**(4), 2**(5)]
#ticks_x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0]
#ticks_y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0]

plt.xticks(ticks_x, va='center')
plt.yticks(ticks_y)

plt.xlim(ticks_x[0] / 2.0, ticks_x[len(ticks_x) - 1] * 2.0)
plt.ylim(ticks_y[0] / 2.0, ticks_y[len(ticks_y) - 1] * 2.0)

ax.xaxis.set_major_formatter(ticker.FuncFormatter(myticks))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(myticks))

ax.tick_params(axis='both', which='major', labelsize=8)

# Adding boundaries from 3 (a)
add_boundaries(ax, 2**(-6), 2**(6), sns.xkcd_rgb['black'])

# Adding intersecition points of boundaries from 3 (a)
add_point(ax, 5.0/32.0, 5.0, '', (0., 0.0), 'grey')
add_point(ax, 20.0/32.0, 20.0, '', (0., 0.0), 'grey')

# Initializing handles for legend
handles = []
legend_labels = ['(b) scalar', '(c) SIMD', '(d) stride=7']

# Adding points from 3 (b)
add_point(ax, 5.0/16.0, 3.0, '(b) $\mathtt{c1}$', (0.1, 0.0), colors[0])
add_point(ax, 5.0/16.0, 2.0, '(b) $\mathtt{c2}$', (0.1, 0.0), colors[0])
h = add_point(ax, 5.0/16.0, 5.0, '(b) $\mathtt{c3}$', (0.1, -0.1), colors[0])
handles.append(h)

# Adding points from 3 (c)
add_point(ax, 5.0/16.0, 10.0, '(c) $\mathtt{c1}$ & $\mathtt{c3}$', (0.1, 0.0), colors[1])
add_point(ax, 5.0/16.0, 8.0, '(c) $\mathtt{c2}$', (0.1, 0.0), colors[1])
h = add_point(ax, 5.0/16.0, 10.0, '', (0.0, 0.0), colors[1])
handles.append(h)

# Adding points from 3 (d)
add_point(ax, 5.0/64.0, 2.5, '(d) $\mathtt{c1}$ & $\mathtt{c3}$', (0.1, 0.0), colors[2])
add_point(ax, 5.0/64.0, 2.0, '(d) $\mathtt{c2}$', (0.1, 0.0), colors[2])
h = add_point(ax, 5.0/64.0, 2.5, '', (0.0, 0.0), colors[2])
handles.append(h)

ax.legend(handles, legend_labels, loc='lower right', fontsize=10)

plt.savefig('data/roofline.png', dpi=300)

