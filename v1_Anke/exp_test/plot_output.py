import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import pandas as pd
# sys.path.append('../../..')
# sys.path.append('../..')
# sys.path.append('..')
from components.plotter_helper import plot_1x3
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


plt.rcParams.update({'font.size':30})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}

if __name__ == "__main__":
    plot_1x3(['3'], ['1'],['0','g','-'],[20],[0,1,2], save='figures/exp3_1x3', labels=['monopolar', 'grounded return', 'active return'])    