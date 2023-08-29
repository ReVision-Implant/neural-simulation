import os
import sys
sys.path.append('..')
sys.path.append('../..')
from bmtk.simulator import bionet
from optparse import OptionParser, BadOptionError, AmbiguousOptionError
from bmtk.simulator.bionet.pyfunction_cache import synaptic_weight
import numpy as np

from mpi4py import MPI
comm = MPI.COMM_WORLD
MPI_RANK = comm.Get_rank()

DEFAULT_LGN = 'inputs/full3_production_3.0sec_SF0.04_TF2.0_ori90.0_c80.0_gs0.5_spikes.trial_0.h5'
DEFAULT_BKG = 'inputs/bkg_spikes_n1_fr1000_dt0.25_100trials.trial_20.h5'


@synaptic_weight
def DirectionRule_others(edge_props, src_node, trg_node):
    nsyn = 1  # edge_props['nsyns']
    sigma = edge_props['weight_sigma']
    src_tuning = src_node['tuning_angle']
    tar_tuning = trg_node['tuning_angle']

    delta_tuning_180 = abs(abs((abs(tar_tuning - src_tuning) % 360.0) - 180.0) - 180.0)
    w_multiplier_180 = np.exp(-(delta_tuning_180 / sigma) ** 2)
    return w_multiplier_180 * edge_props['syn_weight']


@synaptic_weight
def DirectionRule_EE(edge_props, src_node, trg_node):
    nsyn = 1  # edge_props['nsyns']
    sigma = edge_props['weight_sigma']

    src_tuning = src_node['tuning_angle']
    x_src = src_node['x']
    z_src = src_node['z']

    tar_tuning = trg_node['tuning_angle']
    x_tar = trg_node['x']
    z_tar = trg_node['z']

    delta_tuning_180 = abs(abs((abs(tar_tuning - src_tuning) % 360.0) - 180.0) - 180.0)
    w_multiplier_180 = np.exp(-(delta_tuning_180 / sigma) ** 2)

    delta_x = (x_tar - x_src) * 0.07
    delta_z = (z_tar - z_src) * 0.04

    theta_pref = tar_tuning * (np.pi / 180.)
    xz = delta_x * np.cos(theta_pref) + delta_z * np.sin(theta_pref)
    sigma_phase = 1.0
    phase_scale_ratio = np.exp(- (xz ** 2 / (2 * sigma_phase ** 2)))

    # To account for the 0.07 vs 0.04 dimensions. This ensures
    # the horizontal neurons are scaled by 5.5/4 (from the midpoint
    # of 4 & 7). Also, ensures the vertical is scaled by 5.5/7. This
    # was a basic linear estimate to get the numbers (y = ax + b).
    theta_tar_scale = abs(abs(abs(180.0 - abs(tar_tuning) % 360.0) - 90.0) - 90.0)
    phase_scale_ratio = phase_scale_ratio * (5.5 / 4 - 11. / 1680 * theta_tar_scale)

    return w_multiplier_180 * phase_scale_ratio * edge_props['syn_weight']


def run(config_file, **opts):
    conf = bionet.Config.from_json(config_file, validate=True, **opts)
    conf.build_env(MPI_RANK)
    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()
    bionet.nrn.quit_execution()


class PassThroughOptionParser(OptionParser):
    def error(self, msg):
        pass

    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                OptionParser._process_args(self, largs, rargs, values)
            except (BadOptionError, AmbiguousOptionError) as e:
                pass


if __name__ == '__main__':
    parser = PassThroughOptionParser()
    parser.add_option('--no-recurrent', dest='use_recurrent', action='store_false', default=True)
    parser.add_option('--no-lgn', dest='use_lgn', action='store_false', default=True)
    parser.add_option('--no-bkg', dest='use_bkg', action='store_false', default=True)
    parser.add_option('--direction-rule', dest='use_dr', action='store_true', default=False)
    parser.add_option('--lgn-file', dest='lgn_file', action='store', type='string', default=DEFAULT_LGN)
    parser.add_option('--bkg-file', dest='bkg_file', action='store', type='string', default=DEFAULT_BKG)
    parser.add_option('--overwrite', dest='overwrite', action='store_true', default=True)
    options, args = parser.parse_args()

    usr_vars = vars(options)

    # format the output folder
    output_name = 'output'
    if not options.use_recurrent:
        output_name += '_norecurrent'
    if not options.use_lgn:
        output_name += '_nolgn'
    if not options.use_bkg:
        output_name += '_nobkg'
    if options.use_dr:
        output_name += '_directionrule'

    '''
    if not options.overwrite and os.path.exists(output_name):
        for i in range(1, 1000):
            new_name = '{}.{:03d}'.format(output_name, i)
            if not os.path.exists(new_name):
                output_name = new_name
                break

    comm.Barrier()    
    '''
    usr_vars['output_name'] = output_name

    usr_vars['rule'] = ''
    if options.use_dr:
        usr_vars['rule'] = '.direction_rule'

    # Needed for when calling script with nrniv -python run_bionet.py ...
    for arg in args:
        if arg.endswith(__file__):
            args.remove(arg)

    config_file = 'config.json' if len(args) == 0 else args[0]
    run(config_file, **usr_vars)