import pandas as pd 
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.helper import plot_2x2, plot_app
import matplotlib


plt.rcParams.update({'font.size':30})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
# plt.rcParams.update({'axes.labelweight':'bold'})
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']
hfont = {'fontname':'Verdana'}
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}

def setup(save):

    fig = plt.figure(figsize=(6,6))
    ax = plt.gca()
    size=400
    image = np.zeros((size,size))
    radius=size/2
    for i in range(len(image)):
        for j in range(len(image)):
            if (i-size/2+.5)**2 + (j-size/2+.5)**2 >= (size/2)**2:
                image[i,j] = float('NaN')

    masked_array = np.ma.masked_where(np.isnan(image), image)

    cdict = {
        'red': (
            (0.0,  0.0, 0.0),
            (1.0,  0.0, 0.0),
        ),
        'green': (
            (0.0,  0.0, 0.0),
            (1.0,  1.0, 1.0),
        ),
        'blue': (
            (0.0,  0.0, 0.0),
            (1.0,  0.0, 0.0),
        )
    }
    cmap = matplotlib.colors.LinearSegmentedColormap('green', cdict)
    cmap.set_bad(color='white')
    im = ax.imshow(masked_array, cmap=cmap, extent=[-radius,radius,-radius,radius], origin = 'lower', vmin=0, vmax=1)
    ax.set_xticks([-91,91])
    ax.set_yticks([-91,91])
    ax.set_title('opposite orientation', y=1.1)
    for x in [-91,91]:
        sz=300
        for y in [-91,91]:
            if x == y == -91:
                ax.scatter(x,y, marker = 'X', s=sz, c='red', edgecolors='black', label = 'stimulating\nelectrode')
            elif x == 91 and y == -91:
                None
            else:
                if x == y == 91:
                    ax.scatter(x,y, marker = 'X', s=sz, c='white', edgecolors='black', label = 'return\nelectrode')
                else:
                    ax.scatter(x,y, marker = 'X', s=sz, c='white', edgecolors='black')

    ax.scatter(-91,91-50, marker = '$1$', s=sz, c='white', edgecolors='white')
    ax.scatter(91,91-50, marker = '$2$', s=sz, c='white', edgecolors='white')

    # handles, labels = ax.get_legend_handles_labels()
    # fig.legend(handles, labels, loc='center', ncol=1, bbox_to_anchor=(1.2,0.5), framealpha=0)
    plt.tight_layout()
    plt.savefig(save+'.svg', bbox_inches='tight')
    plt.savefig(save+'.png', bbox_inches='tight')
    # plt.show()


if __name__ == "__main__":
    # setup(save='figures/exp2_setup')
    # plot_2x2(['3','2'], ['-'], 20, save='figures/exp2_2x2')
    # plot_2x2(['3-','2-'], ['-'], 20,  save='figures/exp2-_2x2')

    plot_app('2', save='figures/exp2_app1')
    plot_app('2-', save='figures/exp2_app2')

