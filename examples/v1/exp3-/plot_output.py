import pandas as pd 
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.helper import plot_2x3

plt.rcParams.update({'font.size':20})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
# plt.rcParams.update({'axes.labelweight':'bold'})
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']
hfont = {'fontname':'Verdana'}
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}

if __name__ == "__main__":
    plot_2x3(['3','3-'], ['-'], save='figures/exp3-_2x3')
    # plot_2x3(['3','3-'], ['-'], save='figures/exp3-_2x3')

