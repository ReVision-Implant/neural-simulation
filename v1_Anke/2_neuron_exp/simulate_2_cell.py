from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_node_internode/axon_3/waveform_0/amplitude_10/simulation_0',
    config_file='config.json',
    network_dir='networks/network_B',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)