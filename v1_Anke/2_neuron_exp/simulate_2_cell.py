from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_basic_both_far_z',
    config_file='config.json',
    network_dir='networks/both_far_z',
    tstop=100.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)