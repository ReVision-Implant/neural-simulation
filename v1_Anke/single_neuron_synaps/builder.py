import numpy as np
from bmtk.builder.networks import NetworkBuilder

cortex = NetworkBuilder('mcortex')
cortex.add_nodes(
    cell_name='Scnn1a_473845048',
    potental='exc',
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params='472363762_fit.json',
    morphology='Scnn1a_473845048_m.swc'
)

cortex.build()
cortex.save_nodes(output_dir='networks/network_1')

thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(
    N=10,
    pop_name='tON',
    potential='exc',
    model_type='virtual'
)

thalamus.add_edges(
    source={'pop_name': 'tON'}, target=cortex.nodes(),
    connection_rule=5,
    syn_weight=0.001,
    delay=2.0,
    weight_function=None,
    target_sections=['basal', 'apical'],
    distance_range=[0.0, 150.0],
    dynamics_params='AMPA_ExcToExc.json',
    model_template='exp2syn'
)

thalamus.build()
thalamus.save_nodes(output_dir='networks/network_1')
thalamus.save_edges(output_dir='networks/network_1')

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

psg = PoissonSpikeGenerator(population='mthalamus')
psg.add(
    node_ids=range(10),  # Have 10 nodes to match mthalamus
    firing_rate=10.0,    # 10 Hz, we can also pass in a nonhomoegenous function/array
    times=(0.0, 3.0)    # Firing starts at 0 s up to 3 s
)
psg.to_sonata('inputs/mthalamus_spikes.h5')