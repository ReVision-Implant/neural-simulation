from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_mid_neuron/waveform_5ms/sim_basic/amplitude_20/90d_x_close',
    config_file='config.json',
    network_dir='networks/mid_neuron_90d_x',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)