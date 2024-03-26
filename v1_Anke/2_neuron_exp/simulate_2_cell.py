from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_waveform_5ms_pause/sim_basic/amplitude_10/sim_close_x',
    config_file='config.json',
    network_dir='networks/1_close_x',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)