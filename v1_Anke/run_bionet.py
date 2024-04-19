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

def fix_axon_peri_multiple_stubs(hobj, num_stubs, stub_lengths, stub_diameters):
    """
    Replace reconstructed axon with multiple stubs. 
    :param hobj: hoc object
    :param num_stubs: Number of stubs to create
    :param stub_lengths: list of lenghts for the stubs
    :param stub_diameters 
    """

    #Delete existing axon sections
    for sec in hobj.axon:
        h.delete_section(sec=sec)

    #Create new axon structure with specified number of stubs
    h.execute(f'create axon[{num_stubs}]', hobj)

    #Set properties for each stub
    for i in range (num_stubs):
        hobj.axon[i].L=stub_lengths[i]
        hobj.axon[i].diam=stub_diameters[i]
        hobj.axonal.append(sec=hobj.axon[i])
        hobj.all.append(sec=hobj.axon[i])

    #Connect axon sections
    for i in range (num_stubs-1):
        hobj.axon[i+1].connect(hobj.axon[i],1,0) #first parameter 1 = where section is connected; value 1 means the connection is made at the end of the section; which is the distal end of the section being connected to, second parameter: connection is made at the proximal end of the section being connected to

    #connect the first stub to the soma
    hobj.axon[0].connect(hobj.soma[0], 0.5, 0) #connect to the middle of soma (1) and beginning axon stub (0)
    #Define the shape of the axon
    h.define_shape()

    return hobj


def set_params_peri_simpl_hh(hobj, biophys_params):
    """
    :param hobj: NEURON's cell object
    :param biophys_params: name of json file with biophys params for cell's model which determine spiking behavior
    :return:
    """
    passive = biophys_params['passive'][0]
    conditions = biophys_params['conditions'][0]
    genome = biophys_params['genome']

    # Set passive properties
    cm_dict = dict([(c['section'], c['cm']) for c in passive['cm']])
    for sec in hobj.all:
        if "axon" not in sec.name():
            sec.Ra = passive['ra']
            sec.cm = cm_dict[sec.name().split(".")[1][:4]]
            sec.insert('pas')

            for seg in sec:
                seg.pas.e = passive["e_pas"]
            

    # Insert channels and set parameters
    for p in genome:
        dend_sections = [s for s in hobj.all if s.name().split(".")[1][:4] == "dend"]
        soma_sections = [s for s in hobj.all if s.name().split(".")[1][:4] == "soma"]
        axon_sections = [s for s in hobj.all if s.name().split(".")[1][:4] == "axon"]
        apic_sections = [s for s in hobj.all if s.name().split(".")[1][:4] == "apic"]

        if p["section"] == "dend":
            for dend_sec in dend_sections:
                if p["mechanism"] != "":
                    dend_sec.insert(p["mechanism"])     
                setattr(dend_sec, p["name"], p["value"])

        elif p["section"] == "apic":
            for apic_sec in apic_sections:
                if p["mechanism"] != "":
                    apic_sec.insert(p["mechanism"])
                setattr(apic_sec, p["name"], p["value"])    
                        

        elif p["section"] == "soma":
            for soma_sec in soma_sections:
                if p["mechanism"] != "":
                    soma_sec.insert(p["mechanism"])
                setattr(soma_sec, p["name"], p["value"])

        
        elif p["section"] == "axon":
            n=0
            for axon_sec in axon_sections:
                if n % 2 != 0:
                    axon_sec.insert("mammalian_spike_Anke")
                    setattr(axon_sec,"gnatbar_mammalian_spike_Anke", 1.5)
                    setattr(axon_sec,"gnapbar_mammalian_spike_Anke", 0.002)
                    setattr(axon_sec, "gkbar_mammalian_spike_Anke", 1.6)

                    axon_sec.Ra = 150
                    axon_sec.cm = 1.0
                    axon_sec.insert("pas")
                    setattr(axon_sec, "g_pas", 0.04)
                    setattr(axon_sec, "e_pas", -70)
                    axon_sec.ena = 55.0
                    axon_sec.ek = -77.0
                    n+=1

                else: 
                    axon_sec.Ra = 150
                    axon_sec.cm = 0.005
                    axon_sec.insert("pas")
                    setattr(axon_sec, "g_pas", 0.0000015)
                    setattr(axon_sec, "e_pas", -70)
                    n+=1

        else:
            io.log_error(f'another section that was not taken into account!! -> check')    

    # Set reversal potentials
    for erev in conditions['erev']:
        soma_sections=[s for s in hobj.all if s.name().split(".")[1][:4] == "soma"]
        axon_sections=[s for s in hobj.all if s.name().split(".")[1][:4] == "axon"]

        if erev["section"] == "soma":
            #io.log_info(f'erev potentials for soma and axons')
            for soma_sec in soma_sections:
                soma_sec.ena = erev["ena"]
                soma_sec.ek = erev["ek"] 

def aibs_perisomatic(hobj, cell, dynamics_params):
    if dynamics_params is not None:
        node_id = cell["node_id"]
        cell_type = cell['pop_name']       
        io.log_info(f'Fixing cell #{node_id}, {cell_type}')
    
        fix_axon_peri_multiple_stubs(hobj, 4, [30,30,30,30], [1,1,1,1])  
        set_params_peri_simpl_hh(hobj, dynamics_params)

    return hobj

def run(config_file, **opts):
    add_cell_processor(aibs_perisomatic, overwrite=True)
    conf = bionet.Config.from_json(config_file, validate=True, **opts)
    conf.build_env()
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