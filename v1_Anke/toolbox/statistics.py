### Imports

import pandas as pd 
import numpy as np
import sys
sys.path.append('..')
sys.path.append('../bio_components')
from spikes_helper import get_spikes, get_grid
from file_helper import get_dirs, format_params
from scipy.ndimage import gaussian_filter
from scipy.stats import norm


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

    out = np.sum(np.multiply(template,image))

    template = np.sum(np.square(template))
    image = np.sum(np.square(image))
    out = out / np.sqrt(image * template)
    
    return out


def get_corr_intra(sigma, exp, pattern, amplitude, networks):
    
    # Initialisation
    corrs = []
    names = []
    grids = []
    
    # Iterate over networks: get grid and then get normxcorr2 with the 2D grids of the previous iterations
    for i, network in enumerate(networks):
        grid = get_grid(get_spikes(**get_dirs(exp, pattern, amplitude, network)))
        grids.append(gaussian_filter(grid,sigma,truncate=4))
        for j in range(i):
            network1 = np.int(networks[i])
            network2 = np.int(networks[j])
            corrs.append(normxcorr2(grids[network1],grids[network2], 'valid'))
            names.append(str(network1)+str(network2))

    return corrs, names
        
def get_corr_inter(exp, patterns, amplitudes, networks):
    
    # Initialisation
    # assert isinstance(electrodes, list) + isinstance(amplitudes, list) + isinstance(exp, list) == 1 # Correlation between either different electrodes or different amplitudes
    names = []
    n_spikes = []

    # Make lists for iteration (if necessary)
    exp, pattern, amplitude, mice = format_params(exp, pattern, amplitude, mice)

    corrs = np.zeros((len(exp)*len(patterns)*len(amplitudes), len(exp)*len(patterns)*len(amplitudes)))
    N =  np.zeros((len(exp)*len(patterns)*len(amplitudes), len(exp)*len(patterns)*len(amplitudes)))

    # Iterate over electrodes/amplitudes: get spikes for all networks and get corrcoef with spike lists of previous iterations
    for i, exp in enumerate(exp):
        for j,pattern in enumerate(patterns):
            for k,amplitude in enumerate(amplitudes):
                n_spikes = get_spikes(**format_params(exp, pattern, amplitudes, mice))[1]
                n = i*len(electrodes)*len(patterns)+j*len(patterns)+k
                for u in range(n):
                    corrs[u,n] = np.corrcoef(n_spikes[n],n_spikes[u])[0,1]
                    N[u,n] = np.sum((n_spikes[u]+n_spikes[n])>0)
                    names.append(str(n)+'_'+str(u))
    
    return corrs, names, N

def get_p_values(inter_corrs, intra_corr, N):

    p_values = np.zeros(np.shape(inter_corrs))
    Z_beta = 1.645
    C_intra = 0.5*np.log( (1+intra_corr) / (1-intra_corr) )
    for (x,y), inter_corr in np.ndenumerate(inter_corrs):
        N_temp = N[x,y]
        C_inter = 0.5*np.log( (1+inter_corr) / (1-inter_corr) )
        Z_alpha = np.sqrt(N_temp-3) * (C_intra-C_inter) - Z_beta
        print(Z_alpha)
        p = norm.sf(np.abs(Z_alpha))
        p_values[x,y] = p
    return p_values 

if __name__ == "__main__":
    ### Setup

    exps_list = [['3'],['3','3-'],['2','2-']]       # Nested list where each element is a list of which experiments to compare; N_comparisons = len(exps_list)
    amplitudes = ['10','20','30']                   # Which amplitudes to compare
    # exps_list = [['3']]
    # amplitudes = ['10']
    patterns = [0,1,2]
    networks = ['0','1','2']                        # Which networks to compare
    stim_type = '-'                                 # Which stimulation type to compare

    ### Initialisation

    intra_corrs = np.zeros((len(exps_list),len(amplitudes)))    # Initialise zero array of size [N_comparisons x N_amplitudes]

    ### Iterate over N_comparisons and N_amplitudes

    for x,exps in enumerate(exps_list):
        for y,amplitude in enumerate(amplitudes):
            
            if exps == ['2','2-']:                              # If the comparison is between experiment 2 and -2; then
                electrodes = [1,2]                              # Set electrodes to only include 1 and 2 
            
            intra_corr = []                                     # Empty list
            
            # Iterate over N_experiments and N_electrodes and append intraclass correlation for each combination
            for j,exp in enumerate(exps):                      
                for k,pattern in enumerate(patterns):
                    intra_corr.append(get_corr_intra(10, exp, pattern, amplitude, networks)[0])

            # Average intra_corr and save
            intra_corr = np.mean(intra_corr)                 
            intra_corrs[x,y] = intra_corr

            # Save inter-class correlations and corresponding p_values in export/         
            corrs, names, N = get_corr_inter(exps, pattern, amplitude, networks)
            p_values = get_p_values(corrs, intra_corr+0.065, N)
            p_values = pd.DataFrame(p_values)
            corrs = pd.DataFrame(corrs)
            corrs = corrs.round(3)
            corrs = corrs.replace(0, '')
            corrs.to_csv('export1/inter_corrs_' + str(x) + '_' + amplitude + '.csv')
            # p_values = p_values.round(3) #                                         
            p_values.to_csv('export1/p_values_' + str(x) + '_' + amplitude + '.csv')

    np.savetxt('export1/intra_corrs.txt', intra_corrs) # Save intra-class correlation in export/

    # electrodes = [1,2,3]
    # amplitudes = [10,20,30]
    # ratios = np.zeros((3,3))
    # ratios_ = np.zeros((3,3))
    # elec_pos = [[-91,91],[91,91],[91,-91],[-91,-91]]


    # for x,electrode in enumerate(electrodes):
    #     for y,amplitude in enumerate(amplitudes):
    #         centroid = get_centroid_cov(**get_params('3', electrode, '-', amplitude, ['0','1','2']),v1=True)[0]
    #         ratio = np.linalg.norm(centroid[[0,2]] - elec_pos[x]) /  (np.linalg.norm(centroid[[0,2]] - elec_pos[x]) + np.linalg.norm(centroid[[0,2]] - elec_pos[3]))
    #         ratios[y,x] = ratio

    #         centroid_ = get_centroid_cov(**get_params('3-', electrode, '-', amplitude, ['0','1','2']),v1=True)[0]
    #         ratio_ = np.linalg.norm(centroid_[[0,2]] - elec_pos[x]) /  (np.linalg.norm(centroid_[[0,2]] - elec_pos[x]) + np.linalg.norm(centroid_[[0,2]] - elec_pos[3]))
    #         ratios_[y,x] = ratio_

    # print(ratios, '\n', ratios_)

    # print(ttest_rel(ratios.flatten(), ratios_.flatten()))


    # n_spikes = []
    # stim_type = 'g'

    # n_spikes_temp = np.array([])
    # for network in networks:
    #     n_spikes_temp = np.append(n_spikes_temp, get_spikes(**get_params(0, 0, 'g', 20, network))[1])
        
    # n_spikes.append(n_spikes_temp)

    # n_spikes_temp = np.array([])
    # for network in networks:
    #     n_spikes_temp = np.append(n_spikes_temp, get_spikes(**get_params(3, 1, '-', 20, network))[1])
    # n_spikes.append(n_spikes_temp)
    # N = np.sum((n_spikes[0]+n_spikes[1])>0)

    # inter_corr = np.corrcoef(n_spikes[0],n_spikes[1])[0,1]

    # intra_corrs = [get_corr_intra(10, 3, 1, 'g', 20, [0,1,2])[0], get_corr_intra(10, 0, 0, 'g', 20, [0,1,2])[0]]
    # intra_corr = np.mean(intra_corrs)
    # p = get_p_values(np.array([[inter_corr]]), intra_corr, np.array([[N]]))
    # print(inter_corr, intra_corr, N, p)



    # corrs = get_corr_inter(3, 1, '-', [10,20,30],[0,1,2])
    # print(corrs)