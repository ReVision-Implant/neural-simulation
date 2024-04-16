from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_axon_4_diam_1/network_C/waveform_4_5ms/amplitude_10/simulation_0',
    config_file='config.json',
    network_dir='networks/network_C',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)