import pandas as pd 
import numpy as np
import sys
sys.path.append('..')
sys.path.append('../bio_components')
from bio_components.helper import plot_3x3, fractions, get_spikes, get_corr_inter, get_corr_intra, get_p_values, get_centroid_cov, get_params
from scipy.stats import ttest_rel

exps_list = [['3'],['3','3-'],['2','2-']]
amplitudes = ['10','20','30']
# exps_list = [['3']]
# amplitudes = ['10']
electrodes = [1,2,3]
networks = ['0','1','2']
stim_type = '-'
intra_corrs = np.zeros((len(exps_list),len(amplitudes)))

for x,exps in enumerate(exps_list):
    for y,amplitude in enumerate(amplitudes):
        if exps == ['2','2-']:
            electrodes = [1,2]
        intra_corr = []
        for j,exp in enumerate(exps):
            for k,electrode in enumerate(electrodes):
                intra_corr.append(get_corr_intra(10, exp, electrode, stim_type, amplitude, networks)[0])
                
    
        intra_corr = np.mean(intra_corr)
        intra_corrs[x,y] = intra_corr
        corrs, names, N = get_corr_inter(exps, electrodes, stim_type, amplitude, networks)
        p_values = get_p_values(corrs, intra_corr+0.065, N)
        p_values = pd.DataFrame(p_values)
        corrs = pd.DataFrame(corrs)
        corrs = corrs.round(3)
        corrs = corrs.replace(0, '')
        corrs.to_csv('export1/inter_corrs_' + str(x) + '_' + amplitude + '.csv')
        # p_values = p_values.round(3)
        p_values.to_csv('export1/p_values_' + str(x) + '_' + amplitude + '.csv') 

np.savetxt('export1/intra_corrs.txt', intra_corrs)

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