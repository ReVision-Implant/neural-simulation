import pandas as pd
from biophys_components.hdr5 import HDR5
import matplotlib.pyplot as plt
import numpy as np

def plot_activity_3d(nodes_dir, electrodes_dir, spikes_dir, save_dir):
    node_pos = HDR5(nodes_dir).get_positions()
    n_spikes = np.zeros((np.shape(node_pos)[0]))
    elec_pos = pd.read_csv(electrodes_dir, sep=' ')
    elec_pos = elec_pos[['pos_x', 'pos_x', 'pos_x']].to_numpy()[0]

    spikes = pd.read_csv(spikes_dir, sep=' ')
    labels = ['X [$\mu m$]', 'Y [$\mu m$]', 'Z [$\mu m$]']
    for ind in spikes.index:
        n_spikes[spikes['node_ids'][ind]] += 1

    fig = plt.figure(figsize=(9,12))
    ax = plt.axes(projection="3d")
    
    active = node_pos[n_spikes!=0,:]
    inactive = node_pos[n_spikes==0,:]

    p = ax.scatter(active[:,0], active[:,1], active[:,2], marker='o', s=20, cmap='cool', c=n_spikes[n_spikes!=0], label='activated neuron')
    ax.scatter(inactive[:,0], inactive[:,1], inactive[:,2], marker='o', s=20, c='k', label='non-activated neuron')
    ax.scatter(elec_pos[0], elec_pos[1], elec_pos[2], marker = 'd', s=100, color = 'r', label = 'electrode')
    ax.view_init(elev=5., azim=0)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.legend(loc='center', bbox_to_anchor=(1, 0.5))
    cbar = plt.colorbar(p, label='# spikes [-]', shrink=0.5, location='bottom', pad=-0.1)
    cbar.set_ticks(range(1,int(max(n_spikes))+1))
    plt.savefig(save_dir, bbox_inches='tight', transparent=True)
    plt.show()


def plot_activity_distance(nodes_dir, electrodes_dir, spikes_dirs, save_dir, legend=None):
    node_pos = HDR5(nodes_dir).get_positions()
    elec_pos = pd.read_csv(electrodes_dir, sep=' ')
    elec_pos = elec_pos[['pos_x', 'pos_x', 'pos_x']].to_numpy()[0]

    r = np.zeros(np.size(node_pos, axis=0))
    for i in range(np.size(r)):
        r[i] = np.sqrt(((elec_pos-node_pos[i,:])**2).sum())

    fig = plt.figure()
    for spikes_dir in spikes_dirs:
        spikes = pd.read_csv(spikes_dir, sep=' ')
        n_spikes = np.zeros((np.shape(node_pos)[0]))
        for ind in spikes.index:
            n_spikes[spikes['node_ids'][ind]] += 1
        plt.scatter(r, n_spikes, s=20, marker='o')
    
    plt.xlabel('distance to electrode [$\mu m$]')
    plt.ylabel('# spikes [-]')
    plt.title('Number of spikes per neuron')
    plt.yticks(range(0,int(np.max(n_spikes)+1)))
    if legend is not None:
        plt.legend(legend, loc='best')
    plt.savefig(save_dir, transparent=True, bbox_inches='tight')
    plt.show()

plot_activity_distance(
    nodes_dir = 'network/slice_nodes.h5',
    electrodes_dir = 'biophys_components/stimulations/50.csv',
    spikes_dirs = ['output_xstim50/spikes.csv'],
    save_dir = 'figures/a_d_20',
    legend = ['I = 20 $\mu A$']
)

plot_activity_distance(
    nodes_dir = 'network/slice_nodes.h5',
    electrodes_dir = 'biophys_components/stimulations/50.csv',
    spikes_dirs = ['output_xstim50/spikes.csv', 'output_xstim/spikes.csv'],
    save_dir = 'figures/a_d_20_50',
    legend = ['I = 20 $\mu A$', 'I = 50 $\mu A$']
)

plot_activity_3d(
    nodes_dir = 'network/slice_nodes.h5',
    electrodes_dir = 'biophys_components/stimulations/50.csv',
    spikes_dir = 'output_xstim/spikes.csv',
    save_dir = 'figures/3d_50',
)

plot_activity_3d(
    nodes_dir = 'network/slice_nodes.h5',
    electrodes_dir = 'biophys_components/stimulations/50.csv',
    spikes_dir = 'output_xstim50/spikes.csv',
    save_dir = 'figures/3d_20',
)

