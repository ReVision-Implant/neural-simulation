import pandas as pd 
import sys
import numpy as np
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.plot import plotting_params, change_coord, get_spikes, get_grid
from bio_components.statistical_analysis import manova, get_grid, centroid_cov
from scipy.ndimage import gaussian_filter
from bio_components.hdf5 import HDF5
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal
from scipy.signal import correlate2d, fftconvolve
from statsmodels.multivariate.manova import MANOVA
import statsmodels.api as sm
from statsmodels.formula.api import ols


def manova(nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dir=None, t_stop=100, radius=400):

    df = pd.DataFrame(columns=['x','y','z','i'])

    for i in range(len(nodes_dirs)):
        nodes_dir = nodes_dirs[i]
        spikes_dir = spikes_dirs[i]
        elec = spikes_dir[12]
        spikes_bg_dir = spikes_bg_dirs[i]
        node_pos = HDF5(nodes_dir).get_positions_v1()
        n_spikes = np.zeros(np.shape(node_pos)[0])
        spikes = pd.read_csv(spikes_dir, sep='\s+')
        for ind in spikes.index:
            if spikes['timestamps'][ind] < t_stop:
                n_spikes[spikes['node_ids'][ind]] += 1
        spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
        for ind in spikes_bg.index:
            n_spikes[spikes_bg['node_ids'][ind]] = max(0, n_spikes[spikes_bg['node_ids'][ind]] - 1)
        for j in range(len(n_spikes)):
            for k in range(np.int(n_spikes[j])):
                if np.abs(node_pos[j,0]+75)<radius and np.abs(node_pos[j,2]+91)<radius:
                    df.loc[len(df)] = [node_pos[j,0],node_pos[j,1],node_pos[j,2],elec]

    # mod = ols('x ~ i', data=df).fit()
    # return sm.stats.anova_lm(mod)

    maov = MANOVA.from_formula('x + z ~ i', data=df)
    return maov.mv_test() 


def kw(exp, electrodes, stim_types, amplitudes, networks=None, t_stop=100, radius=400):

    nodes_dirs = []
    spikes_dirs = []
    spikes_bg_dirs =[]
    
    for amplitude in amplitudes:
        for electrode in electrodes:
            for stim_type in stim_types:
                nodes_dirs_networks = []
                spikes_dirs_networks = []
                spikes_bg_dirs_networks = []
                for network in networks:
                    nodes_dirs_networks.append('networks_25/network' + network + '/v1_nodes.h5')
                    spikes_dirs_networks.append('exp3/output/'+ electrode + '/' + stim_type + '/' + amplitude + '/0' + electrode + stim_type + '_' + amplitude + '_' + network + '/spikes.csv')
                    spikes_bg_dirs_networks.append('exp3/output/bkg/bkg_' + network + '/spikes.csv')
                nodes_dirs.append(nodes_dirs_networks)
                spikes_dirs.append(spikes_dirs_networks)
                spikes_bg_dirs.append(spikes_bg_dirs_networks)

    electrodes = '../bio_components/stimulations/elec_loc/' + electrode + '.csv'

    for i in range(len(nodes_dirs)):
        phi1 = []
        for j in range(len(nodes_dirs[i])):
            nodes_dir = nodes_dirs[i][j]
            spikes_dir = spikes_dirs[i][j]
            spikes_bg_dir = spikes_bg_dirs[i][j]
            elec = spikes_dir[12]

            node_pos = HDF5(nodes_dir).get_positions_v1()
            n_spikes = np.zeros(np.shape(node_pos)[0])
            spikes = pd.read_csv(spikes_dir, sep='\s+')
            for ind in spikes.index:
                if spikes['timestamps'][ind] < t_stop:
                    n_spikes[spikes['node_ids'][ind]] += 1
            spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
            for ind in spikes_bg.index:
                n_spikes[spikes_bg['node_ids'][ind]] = max(0, n_spikes[spikes_bg['node_ids'][ind]] - 1)
            for k in range(len(n_spikes)):
                for l in range(np.int(n_spikes[k])):
                    if np.abs(node_pos[k,0]+75)<radius and np.abs(node_pos[k,2]+91)<radius:
                        phi = np.arctan2(node_pos[k,0],node_pos[k,2])
                        phi1.append(phi)

        for i2 in range(i+1, len(nodes_dirs)):
            phi2 = []
            for j2 in range(len(nodes_dirs[i2])):

                nodes_dir = nodes_dirs[i2][j2]
                spikes_dir = spikes_dirs[i2][j2]
                elec = spikes_dir[12]
                spikes_bg_dir = spikes_bg_dirs[i2][j2]
                node_pos = HDF5(nodes_dir).get_positions_v1()
                n_spikes = np.zeros(np.shape(node_pos)[0])
                spikes = pd.read_csv(spikes_dir, sep='\s+')
                for ind in spikes.index:
                    if spikes['timestamps'][ind] < t_stop:
                        n_spikes[spikes['node_ids'][ind]] += 1
                spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
                for ind in spikes_bg.index:
                    n_spikes[spikes_bg['node_ids'][ind]] = max(0, n_spikes[spikes_bg['node_ids'][ind]] - 1)
                for k2 in range(len(n_spikes)):
                    for l2 in range(np.int(n_spikes[k2])):
                        if np.abs(node_pos[k2,0]+75)<radius and np.abs(node_pos[k2,2]+91)<radius:
                            phi = np.arctan2(node_pos[k2,0],node_pos[k2,2])
                            phi2.append(phi)
            print(spikes_dirs[i][0], spikes_dirs[i2][0], kruskal(phi1,phi2))
    return
    # mod = ols('x ~ i', data=df).fit()
    # return sm.stats.anova_lm(mod)

def manova_centroids(nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dir=None, t_stop=100, radius=400):

    df = pd.DataFrame(columns=['x','y','z','i'])

    for i in range(len(nodes_dirs)):
        nodes_dir = nodes_dirs[i]
        spikes_dir = spikes_dirs[i]
        elec = spikes_dir[12]
        spikes_bg_dir = spikes_bg_dirs[i]
        centroid, cov= centroid_cov(nodes_dir, spikes_dir, spikes_bg_dir, v1=True, radius=400)
        df.loc[len(df)] = [centroid[0],centroid[1],centroid[2],elec]
    print(df)
    df2 = df.loc[df['i']=='2']
    df2 = df2[['x','y','z']] - df2[['x','y','z']].mean()
    print(df2)
    print(np.cov(df2))
    print([np.sqrt(np.cov(df2[['x','y','z']])[i,i]) for i in range(3)])
    maov = MANOVA.from_formula('x + y + z ~ i', data=df)
    return maov.mv_test() 



def normxcorr2(template, image, mode="full"):
    """
    Input arrays should be floating point numbers.
    :param template: N-D array, of template or filter you are using for cross-correlation.
    Must be less or equal dimensions to image.
    Length of each dimension must be less than length of image.
    :param image: N-D array
    :param mode: Options, "full", "valid", "same"
    full (Default): The output of fftconvolve is the full discrete linear convolution of the inputs. 
    Output size will be image size + 1/2 template size in each dimension.
    valid: The output consists only of those elements that do not rely on the zero-padding.
    same: The output is the same size as image, centered with respect to the ‘full’ output.
    :return: N-D array of same dimensions as image. Size depends on mode parameter.
    """

    # If this happens, it is probably a mistake
    if np.ndim(template) > np.ndim(image) or \
            len([i for i in range(np.ndim(template)) if template.shape[i] > image.shape[i]]) > 0:
        print("normxcorr2: TEMPLATE larger than IMG. Arguments may be swapped.")

    template = template - np.mean(template)
    image = image - np.mean(image)

    a1 = np.ones(template.shape)
    # Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
    ar = np.flipud(np.fliplr(template))
    out = fftconvolve(image, ar.conj(), mode=mode)
    
    image = fftconvolve(np.square(image), a1, mode=mode) - \
            np.square(fftconvolve(image, a1, mode=mode)) / (np.prod(template.shape))

    # Remove small machine precision errors after subtraction
    image[np.where(image < 0)] = 0

    template = np.sum(np.square(template))
    with np.errstate(divide='ignore',invalid='ignore'): 
        out = out / np.sqrt(image * template)

    # Remove any divisions by 0 or very close to 0
    out[np.where(np.logical_not(np.isfinite(out)))] = 0
    
    return out

def corr_anal():
    grid1 = get_grid(**plotting_params('3', '1', '-', '30', ['0']),size=400)
    grid2 = get_grid(**plotting_params('3', '2', '-', '30', ['0']),size=400)
    grid3 = get_grid(**plotting_params('3', '3', '-', '30', ['0']),size=400)
    grid4 = get_grid(**plotting_params('3', '1', '-', '30', ['2']),size=400)

    sigma = 1
    grid1[:,:,1] = gaussian_filter(grid1[:,:,1],sigma,truncate=4) / np.max(gaussian_filter(grid1[:,:,1],sigma,truncate=4))
    grid2[:,:,1] = gaussian_filter(grid2[:,:,1],sigma,truncate=4) / np.max(gaussian_filter(grid2[:,:,1],sigma,truncate=4))
    grid3[:,:,1] = gaussian_filter(grid3[:,:,1],sigma,truncate=4) / np.max(gaussian_filter(grid3[:,:,1],sigma,truncate=4))
    grid4[:,:,1] = gaussian_filter(grid4[:,:,1],sigma,truncate=4) / np.max(gaussian_filter(grid4[:,:,1],sigma,truncate=4))
    fig, axs = plt.subplots(1,4, sharex=True, sharey=True, figsize=(20,9))
    axs[0].imshow(grid1,extent=[0,400,400,0], origin = 'lower')
    axs[1].imshow(grid2,extent=[0,400,400,0], origin = 'lower')
    axs[2].imshow(grid3,extent=[0,400,400,0], origin = 'lower')
    grid_plot = np.zeros((400,400,3))
    grid_plot[:,:,0] = grid1[:,:,1]
    grid_plot[:,:,1] = grid2[:,:,1]
    grid_plot[:,:,2] = grid3[:,:,1]
    axs[-1].imshow(grid_plot,extent=[0,400,400,0], origin = 'lower')
    fig.tight_layout()

    # grid1[:,:,1] -= np.average(grid1[:,:,1])
    # grid2[:,:,1] -= np.average(grid2[:,:,1])
    # grid1 = grid1[:,:,1]
    # grid2 = grid2[:,:,1]
    # grid2[grid2<0]=0
    # grid2[grid2>0]=1
    # plt.figure()
    # plt.imshow(grid2)
    # print(np.sum(grid1[:,:,1]*grid2[:,:,1]) / np.sqrt( np.sum(grid1[:,:,1]**2) * np.sum(grid2[:,:,1]**2) ) )
    print(normxcorr2(grid1[:,:,1],grid2[:,:,1], mode='valid'))
    print(normxcorr2(grid1[:,:,1],grid3[:,:,1], mode='valid'))
    print(normxcorr2(grid2[:,:,1],grid3[:,:,1], mode='valid'))
    print(normxcorr2(grid1[:,:,1],grid4[:,:,1], mode='valid'))
    plt.show()

if __name__ == "__main__":
    res=0

    # res = do_manova(
    #     amplitudes=['20'],
    #     electrodes=['2','3'],
    #     stim_types=['-'],
    #     networks=['0','1','2'],
    #     radius=91
    # )

    # kw('3', ['1','2','3'], ['-'], ['20'], ['0','1','2'], radius=91)
    
    # res = manova_centroids(**plotting_params('3', ['1', '2'], '-', '30', ['0','1','2']))
    
    # print(res)
    corr_anal()