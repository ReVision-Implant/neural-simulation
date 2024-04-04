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
pos_neuron_1=[-11, 60, 12.5]
pos_neuron_2=[-15, 60, 12.5]

net.add_nodes(
    N =2,
    pop_name='e23Cux2',
    positions=np.row_stack((pos_neuron_1, pos_neuron_2)),
    #rotation_angle_xaxis = (np.pi)/2,
    #rotation_angle_yaxis = 0,
    #rotation_angle_zaxis = 0,
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params="487661754_fit.json",
    morphology="Cux2-CreERT2_Ai14-211772.05.02.01_496085150_m.swc"
)

net.add_edges(
    source={'pop_name': 'e23Cux2'}, target={'pop_name': 'e23Cux2'},
    connection_rule=connector,
    syn_weight=100,
    target_sections=['basal', 'apical'],
    distance_range=[0.0,150.0], # determines where on the post-syn cell to place the synapse. The placement is random within the given section and range
    delay=1.6,
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

net.build()
net.save_nodes(output_dir='networks/network_1.5')
net.save_edges(output_dir='networks/network_1.5')
