#This script builds a small population of 100 neurons. The goal is to test different types of axons and parameters

#Option 1: changing the direction of the axon with aibs_perisomatic_directed

import numpy as np

from bmtk.builder.networks import NetworkBuilder
from bmtk.builder.auxi.node_params import positions_columinar, xiter_random
from bmtk.builder.auxi.edge_connectors import distance_connector

net = NetworkBuilder("small_network")
net.add_nodes(
    N=100,
    pop_name='Scnn1a',
    positions=positions_columinar(N=100, center=[0, 50.0, 0], max_radius=30.0, height=100.0),
    rotation_angle_yaxis=xiter_random(N=100, min_x=0.0, max_x=2*np.pi),
    rotation_angle_zaxis=xiter_random(N=100, min_x=0.0, max_x=2*np.pi),
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    #model_processing='aibs_perisomatic_directed', #will find a general direction of the extisting axon before it is cut and use that direction adding the stub
    dynamics_params='472363762_fit.json',
    morphology='Scnn1a_473845048_m.swc'
)


net.add_edges(
    source={'pop_name': 'Scnn1a'}, target={'pop_name': 'Scnn1a'},
    connection_rule=distance_connector,
    connection_params={'d_weight_min': 0.0, 'd_weight_max': 0.34, 'd_max': 50.0, 'nsyn_min': 0, 'nsyn_max': 10},
    syn_weight=2.0e-04,
    distance_range=[30.0, 150.0],
    target_sections=['basal', 'apical', 'soma'],
    delay=2.0,
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')
#net.save_nodes(output_dir='v1_Anke/small_population_exp/network_directed_axons')
#net.save_edges(output_dir='v1_Anke/small_population_exp/network_directed_axons')