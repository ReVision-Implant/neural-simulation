import numpy as np
from bmtk.builder.networks import NetworkBuilder

cortex = NetworkBuilder('mcortex')
cortex.add_nodes(
    pop_name='e23Cux2',
    potental='exc',
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params="487661754_fit.json",
    morphology="Cux2-CreERT2_Ai14-211772.05.02.01_496085150_m.swc"
)

cortex.build()
cortex.save_nodes(output_dir='networks/network_2')

thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(
    N=10,
    pop_name='tON',
    potential='exc',
    model_type='virtual'
)

thalamus.add_edges(
    source={'pop_name': 'tON'}, target=cortex.nodes(),
    connection_rule=2,
    syn_weight=0.001,
    delay=2.0,
    weight_function=None,
    target_sections=['basal', 'apical'],
    distance_range=[0.0, 150.0],
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

thalamus.build()
thalamus.save_nodes(output_dir='networks/network_2')
thalamus.save_edges(output_dir='networks/network_2')

