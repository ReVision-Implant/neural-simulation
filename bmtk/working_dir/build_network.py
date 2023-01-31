import numpy as np
from bmtk.builder import NetworkBuilder
from bmtk.builder.auxi.node_params import positions_columinar, xiter_random
from biophys_components.hdr5 import HDR5
from bmtk.builder.auxi.edge_connectors import distance_connector

import numpy.random
numpy.random.seed(10)

n_nodes = 100

slice = NetworkBuilder('slice')
slice.add_nodes(
    N=n_nodes,
    pop_name='Scnn1a',
    positions=positions_columinar(N=n_nodes, center=[0, 0, 0], max_radius=300, height=150.0, plot=True),
    rotation_angle_xaxis=xiter_random(N=n_nodes, min_x=0.0, max_x=2*np.pi),
    rotation_angle_yaxis=xiter_random(N=n_nodes, min_x=0.0, max_x=2*np.pi),
    rotation_angle_zaxis=xiter_random(N=n_nodes, min_x=0.0, max_x=2*np.pi),
    potental='exc',
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params='472363762_fit.json',
    morphology='Scnn1a_473845048_m.swc'
)

slice.add_edges(
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

slice.build()
slice.save_nodes(output_dir='network')
slice.save_edges(output_dir='network')

# f = HDR5('network/slice_nodes.h5')
# print(f.get_positions())
# print(f.get_rotations())