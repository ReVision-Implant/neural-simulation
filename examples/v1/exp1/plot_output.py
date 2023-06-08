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
from bio_components.helper import get_params, get_image, plot_image, plot_3x3, fractions, get_fractions, plot_2x3
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


plt.rcParams.update({'font.size':30})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}

def create_plot1(amplitude, stim_type, save=None):

    ### Initialisation
    images = []
    logmaxs = []
    electrodes = [8,1,2,7,None,3,6,5,4]

    fig, axs = plt.subplots(3,3, sharex=True, sharey=True, figsize=(16,13))
    for i,electrode in enumerate(electrodes):
        if electrode is not None:
            image = get_image(**get_params('1', electrode, stim_type, amplitude, ['0','1','2']))
            images.append(image)
        else:
            images.append(np.zeros((800,800)))

    images = np.array(images)/np.nanmax(images)

    for i,electrode in enumerate(electrodes):
        if electrode is not None:
            im = plot_image(axs.flatten()[i], images[i,:,:], **get_params('1', electrode, stim_type, amplitude, ['0','1','2']))
        else:
            radius = 400
            masked_array = np.ma.masked_where(np.zeros((2*radius,2*radius)) == 0, np.zeros((2*radius,2*radius)))
            cmap = matplotlib.cm.get_cmap("magma").copy()
            cmap.set_bad(color='white')
            axs[1,1].imshow(masked_array, cmap=cmap, extent=[-radius,radius,-radius,radius], origin = 'lower')

    handles, labels = axs[0,0].get_legend_handles_labels()
    axs[1,1].axes.yaxis.set_visible(False)
    axs[1,1].axes.xaxis.set_visible(False)
    # fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, -0.1), framealpha=0)
    axs[1,1].legend(handles, labels, loc='center', ncol=1, framealpha=0)
    fig.tight_layout(h_pad=0.0,w_pad=1.2)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=30, shrink=1)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()

def setup(save):

    fig = plt.figure(figsize=(6,6))
    ax = plt.gca()
    image = np.zeros((800,800))
    size=800
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
    ax.set_xticks([-182,0,182])
    ax.set_yticks([-182,0,182])
    plt.tight_layout()
    for x in [-182,0,182]:
        sz=300
        for y in [-182,0,182]:
            if x == y == 0:
                ax.scatter(x,y, marker = 'X', s=sz, c='red', edgecolors='black', label = 'stimulating\nelectrode')
            else:
                if x == y == 182:
                    ax.scatter(x,y, marker = 'X', s=sz, c='white', edgecolors='black', label = 'return\nelectrode')
                else:
                    ax.scatter(x,y, marker = 'X', s=sz, c='white', edgecolors='black')

    ax.scatter(-182-50,-182-50, marker = '$6$', s=sz, c='white', edgecolors='white')
    ax.scatter(-182-75,0, marker = '$7$', s=sz, c='white', edgecolors='white')
    ax.scatter(-182-50,182+50, marker = '$8$', s=sz, c='white', edgecolors='white')
    ax.scatter(0,182+75, marker = '$1$', s=sz, c='white', edgecolors='white')
    ax.scatter(0,-182-75, marker = '$5$', s=sz, c='white', edgecolors='white')
    ax.scatter(182+50,-182-50, marker = '$4$', s=sz, c='white', edgecolors='white')
    ax.scatter(182+75,0, marker = '$3$', s=sz, c='white', edgecolors='white')
    ax.scatter(182+50,182+50, marker = '$2$', s=sz, c='white', edgecolors='white')

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center', ncol=1, bbox_to_anchor=(1.2,0.5), framealpha=0)
    plt.savefig(save+'.svg', bbox_inches='tight')
    plt.savefig(save+'.png', bbox_inches='tight')
    # plt.show()

if __name__ == "__main__":
    # setup(save='figures/exp1_setup')

    for amplitude in [10,20,30]:
        create_plot1(amplitude, stim_type='-', save='figures/exp1-'+str(amplitude))
        create_plot1(amplitude, stim_type='g', save='figures/exp1g'+str(amplitude))

    # plot_2x3(exps=['1'], stim_types=['g','-'], save='figures/exp1_2x3')

    plot_3x3(exp=1, stim_type='-', save='figures/exp1-_3x3')
    # plot_3x3(exp=1, stim_type='g', save='figures/exp1g')

    # FTA, FNA, AR = get_fractions(1, '-')
    # print(fractions(FTA, FNA, AR, save='figures/exp1_fractions_-'))
    # FTA, FNA, AR = get_fractions(1, 'g')
    # print(fractions(FTA, FNA, AR, save='figures/exp1_fractions_g'))
    