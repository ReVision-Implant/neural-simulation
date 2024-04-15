#make sure config file also includes a csv file
# activate the correct overwrite function if needed !!

from bmtk.simulator import bionet

from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor
from bmtk.simulator.bionet.default_setters.cell_models import fix_axon_peri, set_params_peri
from bmtk.simulator.bionet.io_tools import io
from bmtk.simulator.bionet.nml_reader import NMLTree
from neuron import h
from ipdb import set_trace
import os
import numpy as np
from sklearn.decomposition import PCA


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


def set_params_peri_axon_copy_soma(hobj, biophys_params):
    """Set biophysical parameters for the cell
    :param hobj: NEURON's cell object
    :param biophys_params: name of json file with biophys params for cell's model which determine spiking behavior
    :return:
    """
    passive = biophys_params['passive'][0]
    conditions = biophys_params['conditions'][0]
    genome = biophys_params['genome']
    io.log_info('copy soma to axons!')

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
                #io.log_info(f'test')
        
        elif p["section"] == "axon":
            #io.log_info(f'axon section nothing happens')
            continue

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
            for axon_sec in axon_sections:
                #io.log_info(f'axon_sec', {axon_sec})
                axon_sec.ena = erev["ena"]
                axon_sec.ek = erev["ek"]    

def set_params_peri_active_axon(hobj, biophys_params):
    """Set biophysical parameters for the cell
    :param hobj: NEURON's cell object
    :param biophys_params: name of json file with biophys params for cell's model which determine spiking behavior
    :return:
    """
    passive = biophys_params['passive'][0]
    conditions = biophys_params['conditions'][0]
    genome = biophys_params['genome']
    io.log_info(f'set active axon properties')

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

        elif p["section"] == "axon":               
            for axon_sec in axon_sections:
                if p["mechanism"] == "":
                    setattr(axon_sec, p["name"], p["value"])

                # Insert mechanisms for spiking
                axon_sec.insert("NaTs")  # Transient sodium current
                axon_sec.insert("K_T")    # Transient potassium current

                # Set parameters for spiking mechanisms
                setattr(axon_sec, "gbar_NaTs", 1.0)  
                setattr(axon_sec, "gbar_K_T", 1.0) 

        else:
            io.log_error(f'another section that was not taken into account!! -> check')    

    # Set reversal potentials
    for erev in conditions['erev']:
        sections = [s for s in hobj.all if s.name().split(".")[1][:4] == erev["section"]]
        for sec in sections:
            sec.ena = erev["ena"]
            sec.ek = erev["ek"]

def set_params_peri_5_channel(hobj, biophys_params):
    io.log_info(f'set 5 channel model')

    passive = biophys_params['passive'][0]
    conditions = biophys_params['conditions'][0]
    genome = biophys_params['genome']

    # Set passive properties --> to do: adapt so not set for axon!! -> other values for axon inserted !
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

        elif p["section"] == "axon":               
            for axon_sec in axon_sections:

                # Insert transient Na and K channels
                axon_sec.insert("mammalian_spike") 
                axon_sec.insert("cad")    

                # Set parameters for spiking mechanisms
                axon_sec.Ra = 136.6 #axial resistance in ohm-cm
                axon_sec.cm = 1 #membrane capacitance (µF/cm2)
                setattr(axon_sec, "gnabar_mammalian_spike", 0.420) #maximum sodium conductance (S/cm^2) 
                setattr(axon_sec, "gkbar_mammalian_spike", 0.250) #maximum potassium delayed rectifier conductance
                setattr(axon_sec,"gcabar_mammalian_spike", 0.00075) #maximum calcium conductance
                setattr(axon_sec, "gkcbar_mammalian_spike", 0.00011) #maximum calcium-dependent potassium conductance
                setattr(axon_sec, "depth_cad", 0.01) #calcium pump depth(microns)
                setattr(axon_sec, "taur_cad",1.5) #time constant (msec)

                setattr(axon_sec, "g_pas", 0.0001)   # leakage conductance
                setattr(axon_sec, "e_pas", -65.02)


                setattr(axon_sec, "ena", 61.02) #sodium resting potential (mV)
                setattr(axon_sec, "ek", -101.31) #potassium resting potential (mV)
                


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

def set_params_peri_hh(hobj, biophys_params):
    io.log_info(f'set hh model')

    passive = biophys_params['passive'][0]
    conditions = biophys_params['conditions'][0]
    genome = biophys_params['genome']

    # Set passive properties --> to do: adapt so not set for axon!! -> other values for axon inserted !
    cm_dict = dict([(c['section'], c['cm']) for c in passive['cm']])
    for sec in hobj.all:
        sec.Ra = passive['ra']
        sec.cm = cm_dict[sec.name().split(".")[1][:4]]
        sec.insert('pas')

        if "axon" not in sec.name():  # Check if the section is not an axon
            for seg in sec:
                seg.pas.e = passive["e_pas"]  

        else:
            print('axon!')

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

        elif p["section"] == "axon":               
            for axon_sec in axon_sections:

                # Insert transient Na and K channels
                axon_sec.insert("hh_model") 

                # Set parameters for spiking mechanisms
                setattr(axon_sec, "gnabar_hh_model", 0.04) #maximum sodium conductance (S/cm^2) 
                setattr(axon_sec, "gkbar_hh_model", 0.035) #maximum potassium delayed rectifier conductance
                setattr(axon_sec, "e_pas", -65) # passive membrane potential
                setattr(axon_sec, "g_pas", 0.003)   # leakage conductance
                setattr(axon_sec, "ena", 55) #sodium resting potential (mV)
                setattr(axon_sec, "ek", -65) #potassium resting potential (mV)


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

def get_axon_direction(hobj):
    for sec in hobj.somatic:
        n3d = int(h.n3d(sec=sec))  # get number of n3d points in each section
        soma_end = np.asarray([h.x3d(n3d - 1, sec=sec), h.y3d(n3d - 1, sec=sec), h.z3d(n3d - 1, sec=sec)])
        mid_point = int(n3d / 2)
        soma_mid = np.asarray([h.x3d(mid_point, sec=sec), h.y3d(mid_point, sec=sec), h.z3d(mid_point, sec=sec)])


    for sec in hobj.all:
        section_name = sec.name().split(".")[1][:4]
        if section_name == 'axon':
            n3d = int(h.n3d(sec=sec))  # get number of n3d points in each section
            axon_p3d = np.zeros((n3d, 3))  # to hold locations of 3D morphology for the current section
            for i in range(n3d):
                axon_p3d[i, 0] = h.x3d(i, sec=sec)
                axon_p3d[i, 1] = h.y3d(i, sec=sec)  # shift coordinates such to place soma at the origin.
                axon_p3d[i, 2] = h.z3d(i, sec=sec)

    # Add soma coordinates to the list
    p3d = np.concatenate(([soma_mid], axon_p3d), axis=0)

    # Compute PCA
    pca = PCA(n_components=3)
    pca.fit(p3d)
    unit_v = pca.components_[0]

    mag_v = np.sqrt(pow(unit_v[0], 2) + pow(unit_v[1], 2) + pow(unit_v[2], 2))
    unit_v[0] = unit_v[0] / mag_v
    unit_v[1] = unit_v[1] / mag_v
    unit_v[2] = unit_v[2] / mag_v

    # Find the direction
    axon_end = axon_p3d[-1] - soma_mid
    if np.dot(unit_v, axon_end) < 0:
        unit_v *= -1

    axon_seg_coor = np.zeros((4, 3))
    # unit_v = np.asarray([0,1,0])
    axon_seg_coor[0] = soma_end
    axon_seg_coor[1] = soma_end + (unit_v * 30.) # 30 units away from the soma in the direction of the axon
    axon_seg_coor[2] = soma_end + (unit_v * 60.)
    axon_seg_coor[3] = soma_end + (unit_v * 90.)

    return axon_seg_coor, soma_mid



def aibs_perisomatic(hobj, cell, dynamics_params):
    if dynamics_params is not None:
        node_id = cell["node_id"]
        cell_type = cell['pop_name']       
        io.log_info(f'Fixing cell #{node_id}, {cell_type}')
        
        #fix_axon_peri(hobj)
        fix_axon_peri_multiple_stubs(hobj, 3, [30,30,30], [1,1,1])
        set_params_peri(hobj, dynamics_params)
        #set_params_peri_axon_copy_soma(hobj, dynamics_params)
        #set_params_peri_active_axon(hobj,dynamics_params)
        #set_params_peri_5_channel(hobj, dynamics_params)
        #set_params_peri_hh(hobj, dynamics_params)

        axon_seg_coordin,soma_mid = get_axon_direction(hobj)
        io.log_info(soma_mid)
        io.log_info(axon_seg_coordin)

    return hobj



#here is the code to edit when just running the simulations, above are all the involved functions
add_cell_processor(aibs_perisomatic, overwrite=True)
dir='sim_axon_3_diam_1/network_B/waveform_4_5ms/amplitude_10/simulation_0'

conf=bionet.Config.from_json(dir+'/config.json')
conf.build_env()
net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)
sim.run()

