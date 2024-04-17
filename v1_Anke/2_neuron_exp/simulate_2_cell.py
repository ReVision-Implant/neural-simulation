from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='modfile_Anke/axon_4_diam_1/network_C/no_stim/simulation_0',
    config_file='config.json',
    network_dir='networks/network_C',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)