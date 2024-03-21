import numpy as np
from bmtk.builder.networks import NetworkBuilder

def connector(source, target):
    count=0
    if source['node_id'] == target['node_id']:
        return None
    if source['node_id'] == 1:
        return None
    return 1

net = NetworkBuilder("net")
pos_neuron_1=[7, 101, 44]
pos_neuron_2=[40, 101, 44]

net.add_nodes(
    N =2,
    pop_name='e23Cux2',
    positions=np.row_stack((pos_neuron_1, pos_neuron_2)),
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params="487661754_fit.json",
    morphology="Cux2-CreERT2_Ai14-211772.05.02.01_496085150_m.swc"
)

net.add_edges(
    source={'pop_name': 'e23Cux2'}, target={'pop_name': 'e23Cux2'},
    connection_rule=connector,
    syn_weight=1,
    target_sections=['basal', 'apical'],
    distance_range=[0.0,1.0], # determines where on the post-syn cell to place the synapse. The placement is random within the given section and range
    delay=1.6,
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

net.build()
net.save_nodes(output_dir='2_neuron_net_full_connect')
net.save_edges(output_dir='2_neuron_net_full_connect')
