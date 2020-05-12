#from scipy import special
from pylab import *
# from scipy.integrate import cumtrapz
# from scipy.stats     import linregress
# from scipy.special   import gamma
# import cmap_bwready as colmap
import sys

# modification of plot parameters
rcParams['text.usetex'] = True
rcParams['font.size'] = 38
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = 'Charter'
rcParams['xtick.major.width'] = 2
rcParams['xtick.major.size'] = 8
rcParams['ytick.major.width'] = 2
rcParams['ytick.major.size'] = 8
# rcParams['xtick.labelsize'] = 'small'
# rcParams['xtick.major.pad'] = 40
# rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


# =================================================================================================
# =================================================================================================
# Information
# =================================================================================================

name = str(sys.argv[1])
if(len(sys.argv)>2):
  Realname = str(sys.argv[2])
else:
  Realname = str(sys.argv[1])


# =================================================================================================
# =================================================================================================
# Cosas automaticas
# =================================================================================================

# loads the data ColoradoSpring1_degdis.txt
structural_data = loadtxt('./' + name + '_degdis.txt')
n_path = structural_data[1:-1,1]/np.sum(structural_data[1:-1,1])

# =================================================================================================
# =================================================================================================
# plots the figure
# =================================================================================================

# Creates the figure and axis objects.
fig, ax = subplots(1,1,figsize=(12,8))
fig.subplots_adjust(bottom=0.15, left=0.15, top=0.95, right=0.95)

# Sets the color of the markers in the scatter plot
# colormap = colmap.cmap_bwready(c_max-c_min+1,1)

# Adds a shaded area.
#ax.axvspan(100, 130, facecolor=[0.75,0.75,0.75], edgecolor='None', alpha=0.5)
#ax.axvspan(215, 230, facecolor=[0.75,0.75,0.75], edgecolor='None', alpha=0.5)


# Plots the histogram.
ax.plot(n_path, c='Black', marker='o', ms=15, ls='-', lw=6, mec='None')
ax.plot(n_path, c='Orange', marker='o', ms=12, ls='-', lw=3, mec='None')

# Adds the name of the network.
text(0.95, 0.90,Realname,
     horizontalalignment='right',
     verticalalignment='center',transform=ax.transAxes)
# ax.annotate(name, 0.01, 0.8, fontsize=12)

# # Sets log-log scale
# ax.set_xscale('log')
# ax.set_yscale('log')

# Bounds.
ax.set_xbound(-1,len(n_path)+1)
#ax.set_ybound(10**(floor(log10(min([i for i in p_his if i>0])))),10**(ceil(log10(p_his.max()))))
# if you need space for name
ax.set_ybound(0,0.2)

# Axes label.
ax.set_xlabel(r'Degree',          fontsize=38, labelpad=0)
ax.set_ylabel(r'Degree distribution', fontsize=38, labelpad=0)

# # Legend.
lines, labels = ax.get_legend_handles_labels()
ax.legend(lines, labels, loc='lower left', shadow=False, fancybox=False, prop={'size':38}, frameon=False, numpoints=1)

# Save to file.
fig.savefig('./' + name + '_degree_distribution.pdf')
clf()


