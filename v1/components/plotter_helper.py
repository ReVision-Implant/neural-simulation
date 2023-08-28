import pandas as pd
from hdf5 import HDF5
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import norm


def get_params(exp, electrodes, stim_type, amplitudes, networks):

    # Initialisation
    nodes = []
    spikes = []
    spikes_bg = []
    elec_loc = []

    # Convert strings to lists if necessary
    if not isinstance(electrodes, list):
        electrodes = [electrodes]
    if not isinstance(amplitudes, list):
        amplitudes = [amplitudes]
    if not isinstance(networks, list):
        networks = [networks]

    # Convert integers to strings if necessary
    exp = str(exp) if not isinstance(exp, str) else exp
    electrodes = [str(i) if not isinstance(i, str) else i for i in electrodes]
    amplitudes = [str(i) if not isinstance(i, str) else i for i in amplitudes]
    networks = [str(i) if not isinstance(i, str) else i for i in networks]

    # Networks to use
    density = '25' if exp=='1' else '100'
    radius = 400 if exp=='1' else 200

    # Iterate over all parameters and get paths for each combination
    for amplitude in amplitudes:
        for electrode in electrodes:
            for network in networks:
                nodes.append('networks_' + density + '/network' + network + '/v1_nodes.h5')
                spikes.append('exp' + exp + '/output/'+ electrode + '/' + stim_type + '/' + amplitude + '/0' + electrode + stim_type + '_' + amplitude + '_' + network + '/spikes.csv')
                spikes_bg.append('exp' + exp + '/output/bkg/bkg_' + network + '/spikes.csv')
            elec_loc.append('../bio_components/stimulations/exp' + exp + '/elec_loc/' + electrode + '.csv')

    return dict(nodes_dirs=nodes, spikes_dirs=spikes, spikes_bg_dirs=spikes_bg, electrodes_dirs=elec_loc, radius=radius, v1=True)

class Plotter():

    def __init__(self, nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dirs, radius, v1=None):
        
        self.nodes_dirs = nodes_dirs
        self.spikes_dirs = spikes_dirs
        self.spikes_bg_dirs = spikes_bg_dirs
        self.electrodes_dirs = electrodes_dirs
        self.radius = radius
        self.v1 = v1

def get_fractions(exp, stim_type):
    radius = 400 if exp == 1 else 200
    ### Initialisation
    images = np.zeros((3,8,radius*2,radius*2))
    electrodes = [1,2,3,4,5,6,7,8]
    amplitudes = [10,20,30]
    FTA = np.zeros((3,8))
    FNA = np.zeros((3,8))
    AR = np.zeros((3,8))
    fig, axs = plt.subplots(3,3, sharex=True, sharey=True, figsize=(16,13))
    for column,electrode in enumerate(electrodes):
        for row,amplitude in enumerate(amplitudes):
            images[row,column,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
            node_pos = np.zeros((0,3))
            n_spikes = np.zeros(0)
            for network in ['0','1','2']:
                node_pos_temp, n_spikes_temp = get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))
                node_pos = np.vstack((node_pos, node_pos_temp))
                n_spikes = np.hstack((n_spikes, n_spikes_temp))
            FNA[row,column] = np.sum(n_spikes>0)/np.sum(n_spikes>-1)
            radius = np.sqrt(node_pos[:,0]**2 + node_pos[:,2]**2)
            AR[row, column] = np.average(radius, weights= n_spikes)
    # images = np.array(images)/np.nanmax(images)*0.9
    # images = np.log10(images+0.1)+1

    images = np.array(images)/np.nanmax(images)

    for i in range(np.size(images,0)):
        for j in range(np.size(images,1)):
            image = images[i,j,:,:]/np.nanmax(images)
            FTA[i,j] = np.sum(image>0.0)/np.sum(image!=None)    

    return FTA, FNA, AR

def get_image(radius, nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dirs, depth = None):
    ### Initalisation
    size = int(radius*2)
    grid = np.zeros((size,size))

    ### Get spikes from each simulation and add them in a single grid array
    for i in range(len(nodes_dirs)):
        grid += get_grid(nodes_dirs[i], spikes_dirs[i], spikes_bg_dirs[i], None, radius=radius, depth = depth)

    ### Gaussian smoothing of the image
    sigma = 10
    image = gaussian_filter(grid, sigma, truncate=4) * sigma / len(nodes_dirs)

    ### Set colour outside circle to NaN
    for i in range(len(image)):
        for j in range(len(image)):
            if (i-size/2+.5)**2 + (j-size/2+.5)**2 >= (size/2)**2:
                image[i,j] = float('NaN')

    return image

def plot_image(ax, image, radius, nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dirs, depth=None):

    masked_array = np.ma.masked_where(np.isnan(image), image)

    cdict = {
        'red': (
            (0.0,  0.0, 0.0),
            (1.0,  0.0, 0.0),
        ),
        'green': (
            (0.0,  0.0, 0.0),
            (0.1,  0.3, 0.3),
            (0.2,  0.45, 0.45),
            (0.4,  0.66, 0.66),
            (0.6,  0.8, 0.8),   
            (0.8,  0.9, 0.9),
            (1.0,  1.0, 1.0),
        ),
        'blue': (
            (0.0,  0.0, 0.0),
            (1.0,  0.0, 0.0),
        )
    }
    cmap = matplotlib.colors.LinearSegmentedColormap('green', cdict)
    cmap.set_bad(color='white')

    ### Plotting
    im = ax.imshow(masked_array, cmap=cmap, extent=[-radius,radius,-radius,radius], origin = 'lower', vmin=0, vmax=1)
    if radius == 200:
        ax.set_xticks([-91,91])
        ax.set_yticks([-91,91])
    elif radius == 400:
        ax.set_xticks([-182,0,182])
        ax.set_yticks([-182,0,182])
    ax.tick_params(direction='inout', bottom=True, top=True, left=True, right=True, width=3, length=7)

    ### Electrode positions
    sz=200
    alpha = 0.7
    electrodes_dir = electrodes_dirs[0]
    elec_pos = pd.read_csv(electrodes_dir, sep=' ')
    elec_pos = elec_pos[['pos_x', 'pos_y', 'pos_z']].to_numpy()
    ax.scatter(elec_pos[0,0], elec_pos[0,2], marker = 'X', s=sz, c='red', edgecolors='black', label = 'stimulating\nelectrode', alpha=0.7)
    ax.scatter(elec_pos[1:,0], elec_pos[1:,2], marker = 'X', s=sz, c='white', edgecolors='black', label = 'return\nelectrode', alpha=0.7)

    ### Centroid positions
    centroid = get_centroid_cov(nodes_dirs, spikes_dirs, spikes_bg_dirs, v1=True, radius=radius, depth=depth)[0]
    ax.scatter(centroid[0], centroid[2], marker = 'X', s=sz, c='lime', edgecolors='black', label = 'centroid', alpha=0.7)

    ### Show centroid-to-electrode distance ratio |centroid->electrode0|/(|centroid->electrode0| + |centroid->electrode1|)
    if np.size(elec_pos,0) > 1:
        ratio = np.linalg.norm(centroid[[0,2]] - elec_pos[0,[0,2]]) /  (np.linalg.norm(centroid[[0,2]] - elec_pos[1,[0,2]]) + np.linalg.norm(centroid[[0,2]] - elec_pos[0,[0,2]]))
    # ax.set_xlabel('centroid ratio = ' + str(np.round(ratio,2)))

    return im
    
def get_centroid_cov(nodes_dirs, spikes_dirs, spikes_bg_dirs=None, v1=False, radius=None, electrodes_dirs=None, depth=None):

    # Initialisation
    totals = np.zeros((1,3))
    weights = np.zeros((1,1))

    # Iterate over networks: get spikes+locations in every network
    for i in range(len(nodes_dirs)):
        nodes_dir = nodes_dirs[i]
        spikes_dir = spikes_dirs[i]
        spikes_bg_dir = spikes_bg_dirs[i]

        node_pos, n_spikes = get_spikes(nodes_dir, spikes_dir, spikes_bg_dir, None, radius, v1)

        in_circle = (node_pos[:,0]**2+node_pos[:,2]**2<radius**2)
        if depth is not None:
            in_circle *= (node_pos[:,1] >= depth[0])*(node_pos[:,1] < depth[1])
        n_spikes = n_spikes[in_circle]
        node_pos = node_pos[in_circle,:]
        totals = np.vstack((totals, node_pos))
        weights = np.append(weights, n_spikes)  

    # Get centroid and cov matrix from spikes+locations from previous step
    centroid = np.average(totals, axis=0, weights=weights)
    cov = np.cov(node_pos.T, fweights=n_spikes)
    return centroid, cov, node_pos, n_spikes

def get_grid(nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dirs, radius, v1=True, depth = None):

    if type(nodes_dirs) is list:
        nodes_dirs = nodes_dirs[0]
        spikes_dirs = spikes_dirs[0]
        spikes_bg_dirs = spikes_bg_dirs[0]
    grid = np.zeros((radius*2,radius*2))

    node_pos, n_spikes = get_spikes(nodes_dirs, spikes_dirs, spikes_bg_dirs, radius=radius, v1=True)
     
    circle = np.sqrt(node_pos[:,0]**2 + node_pos[:,2]**2)
    node_pos = np.array(node_pos[circle<radius,:])
    n_spikes = n_spikes[circle<radius]

    if depth is not None:
        layer = (node_pos[:,1] >= depth[0])*(node_pos[:,1] < depth[1])
        node_pos = np.array(node_pos[layer,:])
        n_spikes = n_spikes[layer]

    for node in range(len(node_pos[:,1])):
        grid_el = (np.floor(node_pos[node,[0,1,2]] + radius)).astype(np.int)
        grid[grid_el[2],grid_el[0]] += n_spikes[node]
        # grid[grid_el[2],grid_el[0]] += min(n_spikes[node],5)
    return grid

def get_spikes(nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dirs=None, radius=None, v1=True):
    
    nodes_dir = nodes_dirs[0] if isinstance(nodes_dirs, list) else nodes_dirs
    spikes_dir = spikes_dirs[0] if isinstance(spikes_dirs, list) else spikes_dirs 

    if v1 == True:
        node_pos = HDF5(nodes_dir).get_positions_v1()
    elif v1 == False:
        node_pos = HDF5(nodes_dir).get_positions()

    n_spikes = np.zeros(np.shape(node_pos)[0])

    spikes = pd.read_csv(spikes_dir, sep='\s+')
    for ind in spikes.index:
        if spikes['timestamps'][ind] < 100:
            n_spikes[spikes['node_ids'][ind]] += 1

    if spikes_bg_dirs is not None:
        spikes_bg_dir = spikes_bg_dirs[0] if isinstance(spikes_bg_dirs, list) else spikes_bg_dirs 
        spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
        for ind in spikes_bg.index:
            if spikes_bg['timestamps'][ind] < 100:
                n_spikes[spikes_bg['node_ids'][ind]] = max(0, n_spikes[spikes_bg['node_ids'][ind]] - 1)

    circle = node_pos[:,0]**2 + node_pos[:,2]**2 < radius**2

    return node_pos[circle,:], n_spikes[circle]


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

def get_corr_intra(sigma, exp, electrode, stim_type, amplitude, networks):
    
    # Initialisation
    corrs = []
    names = []
    grids = []
    
    # Iterate over networks: get grid and then get normxcorr2 with the 2D grids of the previous iterations
    for i, network in enumerate(networks):
        grid = get_grid(**get_params(exp, electrode, stim_type, amplitude, network))
        grids.append(gaussian_filter(grid,sigma,truncate=4))
        for j in range(i):
            network1 = np.int(networks[i])
            network2 = np.int(networks[j])
            corrs.append(normxcorr2(grids[network1],grids[network2], 'valid'))
            names.append(str(network1)+str(network2))

    return corrs, names
        
def get_corr_inter(exp, electrodes, stim_type, amplitudes, networks):
    
    # Initialisation
    # assert isinstance(electrodes, list) + isinstance(amplitudes, list) + isinstance(exp, list) == 1 # Correlation between either different electrodes or different amplitudes
    names = []
    n_spikes = []

    # Make lists for iteration (if necessary)
    exp = [exp] if not isinstance(exp, list) else exp
    electrodes = [electrodes] if not isinstance(electrodes, list) else electrodes
    amplitudes = [amplitudes] if not isinstance(amplitudes, list) else amplitudes
    networks = [networks] if not isinstance(networks, list) else networks

    corrs = np.zeros((len(exp)*len(electrodes)*len(amplitudes), len(exp)*len(electrodes)*len(amplitudes)))
    N =  np.zeros((len(exp)*len(electrodes)*len(amplitudes), len(exp)*len(electrodes)*len(amplitudes)))

    # Iterate over electrodes/amplitudes: get spikes for all networks and get corrcoef with spike lists of previous iterations
    for i, exp in enumerate(exp):
        for j,electrode in enumerate(electrodes):
            for k,amplitude in enumerate(amplitudes):
                n_spikes_temp = np.array([])
                for network in networks:
                    n_spikes_temp = np.append(n_spikes_temp, get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))[1])
                n_spikes.append(n_spikes_temp)
                n = i*len(electrodes)*len(amplitudes)+j*len(amplitudes)+k
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
    None