from hdf5 import HDF5
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_activity_3d(nodes_dir, spikes_dir, spikes_bg_dir, electrodes_dir=None, save_dir=None, square_axis=False, v1=False, show=False):
    if v1 == True:
        node_pos = HDF5(nodes_dir).get_positions_v1()
    elif v1 == False:
        node_pos = HDF5(nodes_dir).get_positions()
    n_spikes = np.zeros((np.shape(node_pos)[0]))

    spikes = pd.read_csv(spikes_dir, sep='\s+')
    for ind in spikes.index:
        n_spikes[spikes['node_ids'][ind]] += 1

    if spikes_bg_dir is not None:
        spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
        for ind in spikes_bg.index:
            n_spikes[spikes_bg['node_ids'][ind]] = n_spikes[spikes_bg['node_ids'][ind]] - 1

    fig = plt.figure(figsize=(9,12))
    ax = plt.axes(projection="3d")
    
    active = node_pos[n_spikes!=0,:]
    inactive = node_pos[n_spikes==0,:]

    ax.scatter(active[:,0], active[:,1], active[:,2], marker='o', s=n_spikes[n_spikes!=0], alpha=0.5, cmap='cool', c=n_spikes[n_spikes!=0], label='activated neuron')
    if v1:
        subset = np.random.choice(np.shape(inactive)[0], round(0.03*np.shape(inactive)[0]))
        inactive = inactive[subset,:]
    ax.scatter(inactive[:,0], inactive[:,1], inactive[:,2], marker='o', s=1, c=subset, alpha=.25, label='non-activated neuron')
   
    if electrodes_dir is not None:
        elec_pos = pd.read_csv(electrodes_dir, sep=' ')
        elec_pos = elec_pos[['pos_x', 'pos_y', 'pos_z']].to_numpy()
        ax.scatter(elec_pos[0,0], elec_pos[0,1], elec_pos[0,2], marker = 's', s=50, color = 'r', label = 'electrode')
        ax.scatter(elec_pos[1:,0], elec_pos[1:,1], elec_pos[1:,2], marker = 's', s=50, color = 'k', label = 'electrode')
    
    ax.view_init(elev=5., azim=-90)
    labels = ['X [$\mu m$]', 'Y [$\mu m$]', 'Z [$\mu m$]']
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.legend(loc='center', bbox_to_anchor=(1, 0.5))

    if square_axis:
        # Create cubic bounding box to simulate equal aspect ratio
        X = node_pos[:,0]
        Y = node_pos[:,1]
        Z = node_pos[:,2]
        # max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
        max_range = np.array([X.max(), Y.max(), Z.max()]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')

    if save_dir is not None:
        plt.savefig(save_dir, bbox_inches='tight', transparent=True)
    if show:
        plt.show()

plot_activity_3d('networks_100/network1/v1_nodes.h5', 'exp3/output/1/-/30/01-_30_1/spikes.csv', 'exp3/output/bkg/bkg_1/spikes.csv',v1=True,show=True)


def plot_positions(nodes_dir, save_dir=None):
    node_pos = HDF5(nodes_dir).get_positions_v1()
    labels = ['X [$\mu m$]', 'Y [$\mu m$]', 'Z [$\mu m$]']

    fig = plt.figure(figsize=(9,12))
    ax = plt.axes(projection="3d")

    p = ax.scatter(node_pos[:,0], node_pos[:,1], node_pos[:,2], marker='o', s=5)
    ax.view_init(elev=5., azim=0)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.legend(loc='center', bbox_to_anchor=(1, 0.5))
    if save_dir is not None:
        plt.savefig(save_dir, bbox_inches='tight', transparent=True)
    plt.show()


def plot_activity_distance(nodes_dir, electrodes_dir, spikes_dirs, save_dir=None, legend=None):
    node_pos = HDF5(nodes_dir).get_positions_v1()
    elec_pos = pd.read_csv(electrodes_dir, sep=' ')
    elec_pos = elec_pos[['pos_x', 'pos_y', 'pos_z']].to_numpy()[0]

    r = np.zeros(np.size(node_pos, axis=0))
    for i in range(np.size(r)):
        r[i] = np.sqrt(((elec_pos-node_pos[i,:])**2).sum())

    fig = plt.figure()
    for spikes_dir in spikes_dirs:
        spikes = pd.read_csv(spikes_dir, sep='\s+')
        n_spikes = np.zeros((np.shape(node_pos)[0]))
        for ind in spikes.index:
            n_spikes[spikes['node_ids'][ind]] += 1
        plt.scatter(r[n_spikes!=0], n_spikes[n_spikes!=0], s=20, marker='o')
    
    plt.xlabel('distance to electrode [$\mu m$]')
    plt.ylabel('# spikes [-]')
    plt.title('Number of spikes per neuron')
    plt.yticks(range(0,int(np.max(n_spikes)+1)))
    if legend is not None:
        plt.legend(legend, loc='best')
    if save_dir is not None:
        plt.savefig(save_dir, transparent=True, bbox_inches='tight')
    plt.show()

def plot_activity_2d(nodes_dir, spikes_dir, electrodes_dir=None, save_dir=None, square_axis=True, spikes_bg_dir=None, v1=False, show=False):
    if v1 == True:
        node_pos = HDF5(nodes_dir).get_positions_v1()
    elif v1 == False:
        node_pos = HDF5(nodes_dir).get_positions()
    n_spikes = np.zeros((np.shape(node_pos)[0]))

    spikes = pd.read_csv(spikes_dir, sep='\s+')
    for ind in spikes.index:
        n_spikes[spikes['node_ids'][ind]] += 1

    if spikes_bg_dir is None:
        spikes_bg_dir = 'exp1/output/bkg/bkg_' + spikes_dir[14] + '/spikes.csv',
        
    spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
    for ind in spikes_bg.index:
        n_spikes[spikes_bg['node_ids'][ind]] = n_spikes[spikes_bg['node_ids'][ind]] - 1

    fig = plt.figure(figsize=(9,12))
    ax = fig.add_subplot()

    active = node_pos[n_spikes!=0,:]
    inactive = node_pos[n_spikes==0,:]

    p = plt.scatter(active[:,0], active[:,2], marker='o', s=2*n_spikes[n_spikes!=0], alpha=0.5, cmap='cool', c=n_spikes[n_spikes!=0], label='activated neuron')
    circle = np.sqrt(inactive[:,0]**2 + inactive[:,2]**2)
    inactive = inactive[circle<400,:]
    if v1:
        subset = np.random.choice(np.shape(inactive)[0], round(0.5*np.shape(inactive)[0]))
        inactive = inactive[subset,:]
    plt.scatter(inactive[:,0], inactive[:,2], marker='o', s=1, c='0.05', alpha=0.05, label='non-activated neuron')
    centroid = centroid_cov(nodes_dir, spikes_dir=spikes_dir, spikes_bg_dir=spikes_bg_dir, v1=True)[0]
    plt.scatter(centroid[0], centroid[2], marker = '8', s=50, color = 'b', label = 'centroid')
    cbar = plt.colorbar(p, label='# spikes [-]', orientation='vertical')
    if len(n_spikes[n_spikes!=0]) > 0:
        cbar.set_ticks(range(int(min(n_spikes[n_spikes!=0])),int(max(n_spikes))+1))
    
    if electrodes_dir is not None:
        elec_pos = pd.read_csv(electrodes_dir, sep=' ')
        elec_pos = elec_pos[['pos_x', 'pos_y', 'pos_z']].to_numpy()
        plt.scatter(elec_pos[0,0], elec_pos[0,2], marker = 's', s=100, color = 'r', label = 'electrode')
        plt.scatter(elec_pos[1:,0], elec_pos[1:,2], marker = 's', s=100, color = 'k', label = 'electrode')
    
    labels = ['X [$\mu m$]', 'Y [$\mu m$]', 'Z [$\mu m$]']
    plt.xlabel(labels[0])
    plt.ylabel(labels[2])
    plt.legend(loc='upper right')

    if square_axis:
        plt.axis(xmin=-400, xmax=400, ymin=-400, ymax=400)
        ax.set_aspect('equal', adjustable='box')
        
    if save_dir is not None:
        plt.savefig(save_dir, bbox_inches='tight', transparent=True)
    if show:
        plt.show()


def plot_activity_2d_smooth(nodes_dirs, spikes_dirs, spikes_bg_dirs=None, electrodes_dir=None, save_dir=None, square_axis=True, v1=False, show=False, t_stop=200):
    
    if spikes_dirs[0][3] == '3':
        radius = 200
    if spikes_dirs[0][3] == '1':
        radius = 400
    factor = 1
    size = int(radius*2/factor)
    for i in range(len(nodes_dirs)):
        nodes_dir = nodes_dirs[i]
        if v1 == True:
            node_pos = HDF5(nodes_dir).get_positions_v1()
        elif v1 == False:
            node_pos = HDF5(nodes_dir).get_positions()

        n_spikes = np.zeros((np.shape(node_pos)[0]))

        spikes_dir = spikes_dirs[i]
        spikes_bg_dir = spikes_bg_dirs[i]

        spikes = pd.read_csv(spikes_dir, sep='\s+')
        for ind in spikes.index:
            if spikes['timestamps'][ind] < t_stop:
                n_spikes[spikes['node_ids'][ind]] += 1

        if spikes_bg_dirs is None:
            spikes_bg_dir = 'exp1/output/bkg/bkg_' + spikes_dir[26] + '/spikes.csv'
        if os.path.exists(spikes_bg_dir):
            spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
            for ind in spikes_bg.index:
                n_spikes[spikes_bg['node_ids'][ind]] = max(0, n_spikes[spikes_bg['node_ids'][ind]] - 1)
            
        active = node_pos[n_spikes!=0,:]/factor
        circle = np.sqrt(active[:,0]**2 + active[:,2]**2)
        active = active[circle<size/2,:]
        # active = np.array(active)
        # active = active[active[:,1]>400,:]
        # active = active[active[:,1]<850,:]
        grid = np.zeros((size,size,3))    
        for node in range(len(active[:,1])):
            grid_el = (np.floor(active[node,[0,2]] + size/2)).astype(np.int)
            grid[size-grid_el[1]-1,grid_el[0],1] += n_spikes[n_spikes!=0][node]
    
    grid[:,:,1] = gaussian_filter(grid[:,:,1],5/factor,truncate=4)
    print(np.max(grid[:,:,1], axis=None))
    grid[:,:,1] *= 3
    # grid[:,:,1] /= np.max(grid[:,:,1], axis=None)

    grid_bg = np.zeros((size,size,3))
    for i in range(len(grid)):
        for j in range(len(grid)):
            if (i-size/2)**2 + (j-size/2)**2 > (size/2)**2:
                grid_bg[i,j,:] = 0.99

    grid_bg[:,:,:] = gaussian_filter(grid_bg[:,:,:],1)

    fig = plt.figure(figsize=(9,12))
    ax = fig.add_subplot()
    ax.imshow(grid+grid_bg)
    plt.xticks(np.arange(0,size,size/8), (np.arange(0,size,size/8)-size/2).astype(int))
    plt.yticks(np.arange(0,size,size/8), (np.arange(size,0,-size/8)-size/2).astype(int))

    if electrodes_dir is None:
        electrodes_dir = '../bio_components/stimulations/elec_loc/' + spikes_dir[12] + '.csv'
    elec_pos = pd.read_csv(electrodes_dir, sep=' ')
    elec_pos = change_coord(elec_pos[['pos_x', 'pos_y', 'pos_z']].to_numpy(), factor=factor, size=size)
    ax.scatter(elec_pos[0,0], elec_pos[0,2], marker = '+', s=150, color = 'red', label = 'central electrode')
    ax.scatter(elec_pos[1:,0], elec_pos[1:,2], marker = '+', s=150, color = 'white', label = 'return electrode')

    centroid = change_coord(centroid_cov(nodes_dir, spikes_dir=spikes_dir, spikes_bg_dir=spikes_bg_dir, v1=True)[0],factor=factor, size=size)
    print(centroid)
    ax.scatter(centroid[0], centroid[2], marker = '+', s=150, color = 'b', label = 'centroid')

    plt.title('return electrode: ' + spikes_dir[12] + ', stimulation type: ' + spikes_dir[14] + ', amplitude: ' + spikes_dir[16:18] + ' uA')

    if save_dir is not None:
        if save_dir == True:
            plt.savefig(spikes_dirs[0][0:4]+'/figures/'+spikes_dirs[0][-19:-11]) 

    ax.legend(facecolor='black',framealpha=0.3)
    if show:
        ax.show()

    return grid+grid_bg


def manova(nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dir=None, t_stop=200):

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
                if node_pos[j,0]**2 + node_pos[j,2]**2 < 100**2:
                    df.loc[len(df)] = [node_pos[j,0],node_pos[j,1],node_pos[j,2],elec]

    # mod = ols('x ~ i', data=df).fit()
    # return sm.stats.anova_lm(mod)

    maov = MANOVA.from_formula('x + z ~ i', data=df)
    return maov.mv_test() 

def do_manova(amplitudes, electrodes, stim_types, networks):
    
    nodes = []
    spikes = []
    spikes_bg =[]
    
    for amplitude in amplitudes:
        for electrode in electrodes:
            for stim_type in stim_types:
                for network in networks:

                    nodes.append('networks_25/network' + network + '/v1_nodes.h5')
                    spikes.append('exp1/output/'+ electrode + '/' + stim_type + '/' + amplitude + '/0' + electrode + stim_type + '_' + amplitude + '_' + network + '/spikes.csv')
                    spikes_bg.append('exp1/output/bkg/bkg_' + network + '/spikes.csv')
    
    electrodes = '../bio_components/stimulations/elec_loc/' + electrode + '.csv'

    return manova(nodes_dirs=nodes, spikes_dirs=spikes, spikes_bg_dirs=spikes_bg, electrodes_dir=electrodes)


def subplot(ax, radius, nodes_dirs, spikes_dirs, spikes_bg_dirs, electrodes_dirs, maxlogmax=1):
    
    ### Initalisation
    size = int(radius*2)
    grid = np.zeros((size,size))

    ### Get spikes from each simulation and add them in a single grid array
    for i in range(len(nodes_dirs)):
         grid += get_grid(nodes_dirs[i], spikes_dirs[i], spikes_bg_dirs[i], None, radius=radius)

    ### Gaussian smoothing of the image
    sigma = 10
    grid = gaussian_filter(grid, sigma, truncate=4) * sigma / len(nodes_dirs)
    print(np.min(grid),np.max(grid))
    grid = np.log10(grid+0.1)
    grid += 1
    logmax = np.max(grid)
    grid /= logmax

    ### Set colour outside circle to white
    for i in range(len(grid)):
        for j in range(len(grid)):
            if (i-size/2+.5)**2 + (j-size/2+.5)**2 >= (size/2)**2:
                grid[i,j] = None
    
    masked_array = np.ma.masked_where(grid == None, grid)
    cmap = matplotlib.cm.inferno  # Can be any colormap that you want after the cm

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

    ### Plotting
    im = ax.imshow(masked_array, cmap=cmap, extent=[-radius,radius,-radius,radius], origin = 'lower', vmin=0, vmax=1)
    if radius == 200:
        ax.set_xticks([-91,91])
        ax.set_yticks([-91,91])
    elif radius == 400:
        ax.set_xticks([-182,0,182])
        ax.set_yticks([-182,0,182])

    ### Electrode positions
    sz=100
    electrodes_dir = electrodes_dirs[0]
    elec_pos = pd.read_csv(electrodes_dir, sep=' ')
    elec_pos = elec_pos[['pos_x', 'pos_y', 'pos_z']].to_numpy()
    ax.scatter(elec_pos[0,0], elec_pos[0,2], marker = 'X', s=sz, c='red', edgecolors='black', label = 'central electrode')
    ax.scatter(elec_pos[1:,0], elec_pos[1:,2], marker = 'X', s=sz, c='white', edgecolors='black', label = 'return electrode')

    ### Centroid positions
    centroid = get_centroid_cov(nodes_dirs, spikes_dirs, spikes_bg_dirs, v1=True, radius=radius)[0]
    ax.scatter(centroid[0], centroid[2], marker = 'X', s=sz, c='lime', edgecolors='black', label = 'centroid')

    ### Show centroid-to-electrode distance ratio |centroid->electrode0|/(|centroid->electrode0| + |centroid->electrode1|)
    ratio = np.linalg.norm(centroid[[0,2]] - elec_pos[0,[0,2]]) /  (np.linalg.norm(centroid[[0,2]] - elec_pos[1,[0,2]]) + np.linalg.norm(centroid[[0,2]] - elec_pos[0,[0,2]]))
    # ax.set_xlabel('centroid ratio = ' + str(np.round(ratio,2)))

    return grid, logmax, im