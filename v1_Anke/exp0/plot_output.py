import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import pandas as pd
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.helper import get_params, get_image, plot_image, plot_3x3, fractions, get_fractions, plot_1x3
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


plt.rcParams.update({'font.size':30})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}

if __name__ == "__main__":
    plot_1x3(['3'], ['1'],['0','g','-'],[20],[0,1,2], save='figures/exp3_1x3', labels=['monopolar', 'grounded return', 'active return'])    