from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_waveform_5ms_pause/axon_3_diam_13/amplitude_20/copy_soma_n58',
    config_file='config.json',
    network_dir='networks/neuron_pos_58',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)