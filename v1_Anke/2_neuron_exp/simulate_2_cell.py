from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_internode_tests/network_F/waveform_0/simulation_0',
    config_file='config.json',
    network_dir='networks/network_F',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)