from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_waveform_with_pause_2s/sim_axon_5_diam_1/amplitude_20/conduct_copy_soma',
    config_file='config.json',
    network_dir='networks/far_z',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)