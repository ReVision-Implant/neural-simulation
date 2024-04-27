from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='simulation_9',
    config_file='config.json',
    network_dir='networks/network_1',
    tstop=3000.0, dt=0.1,
    report_vars=['v'],    # Record membrane potential and calcium (default soma)
    spikes_inputs=[('mthalamus', # Name of population which spikes will be generated for
                    'inputs/mthalamus_spikes.h5')],
    include_examples=True,       # Copies components files
    compile_mechanisms=True      # Will try to compile NEURON mechanisms
)