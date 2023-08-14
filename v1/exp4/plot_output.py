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
from bio_components.helper import get_corr_inter, plot_1x2, get_corr_intra

plt.rcParams.update({'font.size':20})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
# plt.rcParams.update({'axes.labelweight':'bold'})
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']
hfont = {'fontname':'Verdana'}
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}

plot_1x2([3,4], [1], ['-'], [20], [0,1,2], labels=['heuristic pulse','rectangular pulse'], save='figures/pulse')
print(get_corr_inter(['3','4'],1,'-',20,['0','1','2']))
print(np.mean([get_corr_intra(10,'3',1,'-',20,['0','1','2'])[0],get_corr_intra(10,'4',1,'-',20,['0','1','2'])[0]]))