from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_mid_neuron/sim_basic/amplitude_20/horizontal_close',
    config_file='config.json',
    network_dir='networks/mid_neuron_horizontal',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)