from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='simulation_basic',
    config_file='config.json',
    network_dir='2_neuron_net_full_connect',
    tstop=100.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)