#This script builds a small population of 100 neurons. The goal is to test different types of axons and parameters

#Option 1: changing the direction of the axon with aibs_perisomatic_directed

module_path='/Users/ankev/Documents/Github/neural-simulation/v1_Anke'
import sys;
sys.path.append(module_path)

import numpy as np

from bmtk.builder.networks import NetworkBuilder
from bmtk.builder.auxi.node_params import positions_columinar, xiter_random
from bmtk.builder.auxi.edge_connectors import distance_connector
from build_files.edge_funcs import connect_cells

net = NetworkBuilder("small_network")
net.add_nodes(
    N=100,
    pop_name='e23Cux2',
    positions=positions_columinar(N=100, center=[0, 50.0, 0], max_radius=30.0, height=100.0),
    rotation_angle_yaxis=xiter_random(N=100, min_x=0.0, max_x=2*np.pi),
    rotation_angle_zaxis=xiter_random(N=100, min_x=0.0, max_x=2*np.pi),
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    #model_processing='aibs_perisomatic',
    #model_processing='aibs_perisomatic_directed', #will find a general direction of the extisting axon before it is cut and use that direction adding the stub
    dynamics_params="487661754_fit.json",
    morphology="Cux2-CreERT2_Ai14-211772.05.02.01_496085150_m.swc"
)


net.add_edges(
    source={'pop_name': 'e23Cux2'}, target={'pop_name': 'e23Cux2'},
    connection_rule=distance_connector,
    connection_params={'d_weight_min': 0.0, 'd_weight_max': 0.34, 'd_max': 50.0, 'nsyn_min': 0, 'nsyn_max': 10},
    syn_weight=0.000182520075015,
    distance_range=[0.0, 200.0],
    target_sections=['basal', 'apical'],
    delay=1.6,
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

net.build()
#net.save_nodes(output_dir='network')
#net.save_edges(output_dir='network')
#net.save_nodes(output_dir='network_directed_axons')
#net.save_edges(output_dir='network_directed_axons')