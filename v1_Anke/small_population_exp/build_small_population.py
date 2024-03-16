#This script builds a small population of 100 neurons. The goal is to test different types of axons and parameters

module_path='/Users/ankev/Documents/Github/neural-simulation/v1_Anke'
import sys;
sys.path.append(module_path)

import numpy as np

from bmtk.builder.networks import NetworkBuilder
from bmtk.builder.auxi.node_params import positions_columnar, xiter_random
from bmtk.builder.auxi.edge_connectors import distance_connector
from build_files.edge_funcs import connect_cells
from build_files.node_funcs import generate_random_positions

import sys;
module_path='/Users/ankev/Documents/Github/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path)
from mask_small import apply_mask_small

net = NetworkBuilder("small_network")
net.add_nodes(
    N=137,
    pop_name='e23Cux2',
    #positions = positions_columnar(N=137, center=[0, 50.0, 0], max_radius=65.0, height=80.0),
    positions_no_mask = generate_random_positions(N, [40,120], [0,65]),
    positions=apply_mask_small(positions_no_mask),
    rotation_angle_yaxis=xiter_random(N=100, min_x=0.0, max_x=2*np.pi),
    rotation_angle_zaxis=xiter_random(N=100, min_x=0.0, max_x=2*np.pi),
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params="487661754_fit.json",
    morphology="Cux2-CreERT2_Ai14-211772.05.02.01_496085150_m.swc"
)


net.add_edges(
    source={'pop_name': 'e23Cux2'}, target={'pop_name': 'e23Cux2'},
    connection_rule=distance_connector, #differs to the bigger model; connections based on probability based on distance
    connection_params={'d_weight_min': 0.0, 'd_weight_max': 0.34, 'd_max': 50.0, 'nsyn_min': 0, 'nsyn_max': 10},
    syn_weight=0.000182520075015,
    distance_range=[0.0, 200.0],
    target_sections=['basal', 'apical'],
    delay=1.6,
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

net.build()
net.save_nodes(output_dir='network_mask_dense')
net.save_edges(output_dir='network_mask_dense')
