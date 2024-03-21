from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='simulation_basic_both_far',
    config_file='config.json',
    network_dir='both_far',
    tstop=100.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)