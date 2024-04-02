from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir='sim_waveform_5ms_pause/axon_2_diam_1/test_n58',
    config_file='config.json',
    network_dir='networks/neuron_pos_58',
    tstop=3000.0, dt=0.025,
    report_vars=['v'], # Record membrane potential
)