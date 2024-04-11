from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_axon_2_diam_1/network_1/waveform_long_pause/simulation_0',
    config_file='config.json',
    network_dir='networks/network_1',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)