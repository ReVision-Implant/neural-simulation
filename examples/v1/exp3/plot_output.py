import pandas as pd 
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.helper import plot_3x3, fractions, plot_1x3, plot_depth
import matplotlib

plt.rcParams.update({'font.size':30})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
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
    if save == 'figures/exp3_setup2': 
        plt.title('original orientation', y=1.1)
    for x in [-91,91]:
        sz=300
        for y in [-91,91]:
            if x == y == -91:
                ax.scatter(x,y, marker = 'X', s=sz, c='red', edgecolors='black', label = 'stimulating\nelectrode')
            else:
                if x == y == 91:
                    ax.scatter(x,y, marker = 'X', s=sz, c='white', edgecolors='black', label = 'return\nelectrode')
                else:
                    ax.scatter(x,y, marker = 'X', s=sz, c='white', edgecolors='black')

    ax.scatter(-91,91+50, marker = '$1$', s=sz, c='white', edgecolors='white')
    ax.scatter(91,91+50, marker = '$2$', s=sz, c='white', edgecolors='white')
    ax.scatter(91,-91+50, marker = '$3$', s=sz, c='white', edgecolors='white')
    plt.tight_layout()
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center', ncol=1, bbox_to_anchor=(1.2,0.5), framealpha=0)
    plt.savefig(save+'.svg', bbox_inches='tight')
    plt.savefig(save+'.png', bbox_inches='tight')
    # plt.show()


if __name__ == "__main__":

    # setup('figures/exp3_setup')
    # setup('figures/exp3_setup2')
    FTA, FNA, AR = plot_3x3(exp=3, stim_type='-', save='figures/exp3_3x3')
    print(fractions(FTA, FNA, AR, save='figures/exp3_fractions'))
    # plot_1x2(['3'],['1'],['g','-'],['20'],['0','1','2'], save='figures/exp3_1x2',labels=['grounded return', 'active return'])
    # plot_depth(3,1,'-',30,[0,1,2], depths=[50,425,850], save='figures/depth', labels=['L1-6','L1-4','L5-6'])