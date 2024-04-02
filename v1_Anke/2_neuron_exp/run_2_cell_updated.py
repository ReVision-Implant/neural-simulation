#make sure config file also includes a csv file
# activate the correct overwrite function if needed !!

from bmtk.simulator import bionet

from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor
from bmtk.simulator.bionet.default_setters.cell_models import fix_axon_peri, set_params_peri
from bmtk.simulator.bionet.io_tools import io
from neuron import h
from bmtk.analyzer.compartment import plot_traces
# from ipdb import set_trace

import matplotlib.pyplot as plt

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
    hobj.axon[0].connect(hobj.soma[0], 1, 0) #0.5 means that the connection is made at the midpoint of the soma

    #Define the shape of the axon
    h.define_shape()

    return hobj


def set_params_peri_axon_copy_soma(hobj, biophys_params):
    """Set biophysical parameters for the cell
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
            #io.log_info(f'dyn param dend')
            for dend_sec in dend_sections:
                if p["mechanism"] != "":
                    dend_sec.insert(p["mechanism"])     
                setattr(dend_sec, p["name"], p["value"])

        elif p["section"] == "apic":
            #io.log_info(f'dyn param apic')
            for apic_sec in apic_sections:
                if p["mechanism"] != "":
                    apic_sec.insert(p["mechanism"])
                setattr(apic_sec, p["name"], p["value"])    
                        

        elif p["section"] == "soma":
            #io.log_info(f'soma & axon section param set')
            for soma_sec in soma_sections:
                if p["mechanism"] != "":
                    soma_sec.insert(p["mechanism"])
                setattr(soma_sec, p["name"], p["value"])        
            for axon_sec in axon_sections:
                if p["mechanism"] != "":
                    axon_sec.insert(p["mechanism"])
                setattr(axon_sec, p["name"], p["value"])
        
        elif p["section"] == "axon":
            #io.log_info(f'axon section nothing happens')
            continue

        else:
            io.log_error(f'another section that was not taken into account!! -> check')    

    # Set reversal potentials
    for erev in conditions['erev']:
        soma_sections=[s for s in hobj.all if s.name().split(".")[1][:4] == "soma"]
        axon_sections=[s for s in hobj.all if s.name().split(".")[1][:4] == "axon"]
        #dend_sections = [s for s in hobj.all if s.name().split(".")[1][:4] == "dend"]
        #apic_sections = [s for s in hobj.all if s.name().split(".")[1][:4] == "apic"]

        #if erev["section"] == "dend":  -> left out because dendrites always passive in v1 model bmtk
        #    for dend_sec in dend_sections:
        #        dend_sec.ena=erev["ena"]
        #        dend_sec.ena=erev["ek"]

        if erev["section"] == "soma":
            #io.log_info(f'erev potentials for soma and axons')
            for soma_sec in soma_sections:
                soma_sec.ena = erev["ena"]
                soma_sec.ek = erev["ek"]
            for axon_sec in axon_sections:
                axon_sec.ena = erev["ena"]
                axon_sec.ek = erev["ek"]    

def aibs_perisomatic(hobj, cell, dynamics_params):
    if dynamics_params is not None:
        node_id = cell["node_id"]
        cell_type = cell['pop_name']       
        io.log_info(f'Fixing cell #{node_id}, {cell_type}')
        
        # fix_axon_peri(hobj)
        fix_axon_peri_multiple_stubs(hobj, 10, [30,30,30,30,30,30,30,30,30,30], [12,12,12,12,12,12,12,12,12,12])
        #set_params_peri(hobj, dynamics_params)
        set_params_peri_axon_copy_soma(hobj, dynamics_params)

    return hobj

def plot_axon_vm(output_dir='output', pop_name='net'):
    import h5py
    import numpy as np
    import matplotlib.pyplot as plt

    plt.figure()
    v_report = h5py.File(f'{output_dir}/v_report.h5', 'r')
    node_ids = v_report[f'report/{pop_name}/mapping/node_ids'][()]
    idx_ptr = v_report[f'report/{pop_name}/mapping/index_pointer'][()]
    seg_types = v_report[f'report/{pop_name}/mapping/seg_types']
    data = v_report[f'report/{pop_name}/data']
    for idx, node_id in enumerate(node_ids):
        cell_beg, cell_end = idx_ptr[idx], idx_ptr[idx+1]
        cell_segs = seg_types[cell_beg:cell_end]
        axon_idxs = np.argwhere(cell_segs == 2).flatten() + cell_beg
        print(f'{node_id} - axon indices {axon_idxs}')
        mean_axon_vms = np.mean(data[:, axon_idxs], axis=1)
        plt.plot(mean_axon_vms, label=node_id)
    
    plt.title('Axon Traces')
    plt.legend()
    # plt.show()


def plot_axon_extracell_input(output_dir='output', pop_name='net'):
    import h5py
    import numpy as np
    

    plt.figure()
    v_report = h5py.File(f'{output_dir}/e_report.h5', 'r')
    node_ids = v_report[f'report/{pop_name}/mapping/node_ids'][()]
    idx_ptr = v_report[f'report/{pop_name}/mapping/index_pointer'][()]
    seg_types = v_report[f'report/{pop_name}/mapping/seg_types']
    data = v_report[f'report/{pop_name}/data']
    for idx, node_id in enumerate(node_ids):
        cell_beg, cell_end = idx_ptr[idx], idx_ptr[idx+1]
        cell_segs = seg_types[cell_beg:cell_end]
        axon_idxs = np.argwhere(cell_segs == 2).flatten() + cell_beg
        print(f'{node_id} - axon indices {axon_idxs}')
        mean_axon_vms = np.mean(data[:, axon_idxs], axis=1)
        plt.plot(mean_axon_vms, label=node_id)
    
    plt.title('Extracellular input (axon)')
    plt.legend()
    
add_cell_processor(aibs_perisomatic, overwrite=True)




dir='sim_waveform_5ms_pause/axon_10_diam_1/amplitude_20/conduct_copy_soma_n58'
output= dir+'/output'

conf=bionet.Config.from_json(dir+'/config.json')
conf.build_env()
net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)
sim.run()

# plot_traces(config_file='config.json', report_name='v_report', population='net', show=False)
plot_axon_extracell_input(output_dir = output)
plot_axon_vm(output_dir=output)
plt.show()
