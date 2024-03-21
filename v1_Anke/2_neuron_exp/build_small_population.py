
import numpy as np
from bmtk.builder.networks import NetworkBuilder


net = NetworkBuilder("2_neuron_network")
pos_neuron_1=[7,101,44]
pos_neuron_2=[40,101,44]
net.add_nodes(
    N=2,
    pop_name='e23Cux2',
    positions=np.row_stack((pos_neuron_1,pos_neuron_2)),
    #rotation_angle_xais=xiter_random(N=137, min_x=0.0, max_x=2*np.pi),
    #rotation_angle_yaxis=xiter_random(N=137, min_x=0.0, max_x=2*np.pi),
    #rotation_angle_zaxis=xiter_random(N=137, min_x=0.0, max_x=2*np.pi),
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params="487661754_fit.json",
    morphology="Cux2-CreERT2_Ai14-211772.05.02.01_496085150_m.swc"
)


net.add_edges(
    source={'pop_name': 'e23Cux2'}, target={'pop_name': 'e23Cux2'},
    connection_rule=1,
    syn_weight=1,
    target_sections=['basal', 'apical'],
    delay=1.6,
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

net.build()
net.save_nodes(output_dir='2_neuron_net_full_connect')
net.save_edges(output_dir='2_neuron_net_full_connect')
