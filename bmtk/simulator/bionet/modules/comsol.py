import os
import math
import pandas as pd
import numpy as np
import six
from neuron import h

from scipy.interpolate import NearestNDInterpolator as NNip
from scipy.interpolate import LinearNDInterpolator as Lip

from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
from bmtk.simulator.bionet.modules.xstim_waveforms import stimx_waveform_factory

class ComsolMod(SimulatorMod):
    """ 
    __init__: COMSOL output .txt file is loaded as pandas dataframe and then used to set up nearest neighbour (NN) interpolation object to create interpolation map
    :param comsol_file: path of .txt file. Coordinates in [um], potentials in [mV], timepoints in [s].
    :param waveform: path of .csv file as created with examples/bio_components/waveform.py
                    If specified, comsol_file should contain output from stationary study and waveform is defined through this parameter.
                    If not specified, comsol_file should contain output from time-dependent study.

    initialise: An interpolation map is defined of every segment and stored in dictionary self._NN, done iteratively for every cell/gid.
    The interpolation map points (the center of) every segment to its NN. It is calculated once here and then used in every step. 
    Next, the COMSOL output is also interpolated in time to match the timesteps in BMTK.

    step: The interpolation map is used to point each segment to its NN and find the corresponding voltage value in the comsol df.
    """

    def __init__(self, comsol_files, waveforms=None, amplitudes=1, 
                 cells=None, set_nrn_mechanisms=True, node_set=None):
        
        if waveforms is None:
            self._comsol_files = comsol_files 
            self._waveforms = waveforms
            self._amplitudes = amplitudes

        else:
            self._comsol_files = comsol_files if type(comsol_files) is list else [comsol_files]
            self._nb_files = len(self._comsol_files) 
            self._waveforms = waveforms if type(waveforms) is list else [waveforms]
            _amplitudes = amplitudes if type(amplitudes) is list else [amplitudes]
            self._amplitudes = _amplitudes*len(self._comsol_files) if len(_amplitudes) == 1 else _amplitudes
            
            try:
                assert self._nb_files == len(self._comsol_files) == len(self._waveforms) == len(self._amplitudes)
            except AssertionError:
                print("AssertionError: comsol_files, waveforms, and amplitudes have a different length.")

            self._data = [None]*self._nb_files

        self._set_nrn_mechanisms = set_nrn_mechanisms
        self._cells = cells
        self._local_gids = []
        

    def initialize(self, sim):
        if self._cells is None:
            # if specific gids not listed just get all biophysically detailed cells on this rank
            self._local_gids = sim.biophysical_gids
        else:
            # get subset of selected gids only on this rank
            self._local_gids = list(set(sim.local_gids) & set(self._all_gids))
        
        if self._waveforms is None:
            self._data =  self.load_comsol(self._comsol_files)
            self._NNip = NNip(self._data[['x','y','z']], np.arange(len(self._data['x'])))
            self._NN = {}
            
            for gid in self._local_gids:
                cell = sim.net.get_cell_gid(gid)
                cell.setup_xstim(self._set_nrn_mechanisms)
                
                r05 = cell.seg_coords.p05               # Get position of segment centre
                self._NN[gid] = self._NNip(r05.T)       # Spatial interpolation
                
            # Temporal interpolation
            timestamps_comsol = np.array(list(self._data)[3:], dtype=float)[:,0]
            timestamps_bmtk = np.arange(timestamps_comsol[0], timestamps_comsol[-1]+sim.dt, sim.dt)
            self._data_temp = np.zeros((self._data.shape[0], len(timestamps_bmtk)))           
            for i in range(self._data.shape[0]):                                 
                self._data_temp[i,:] = np.interp(timestamps_bmtk, timestamps_comsol, self._data.iloc[i,3:]).flatten()                                                          
            self._data = self._data_temp*self._amplitudes
            self._period = int(timestamps_bmtk[-1]/sim.dt)

        else:
            self._Lip = [None]*self._nb_files
            self._L = [None]*self._nb_files
            self._arr = [None]*self._nb_files

            for i in range(self._nb_files):
                self._data[i] =  self.load_comsol(self._comsol_files[i])
                self._waveforms[i] = stimx_waveform_factory(self._waveforms[i])                                                           
                self._arr[i] = self._data[0].to_numpy().flatten()

                self._Lip[i] = Lip(self._data[i][['x','y','z']], self._data[i][0])
                self._L[i] = {}

            for gid in self._local_gids:
                cell = sim.net.get_cell_gid(gid)
                cell.setup_xstim(self._set_nrn_mechanisms)

                r05 = cell.seg_coords.p05               # Get position of middle of segment
                for i in range(self._nb_files):          # Spatial interpolation
                    self._L[i][gid] = self._Lip[i](r05.T)         

                
    def step(self, sim, tstep):
        for gid in self._local_gids:
            cell = sim.net.get_cell_gid(gid)

            if self._waveforms is None:
                NN = self._NN[gid]               # vector that points each node of the cell (BMTK) to the nearest COMSOL node
                tstep = tstep % self._period     # In case of periodic stimulation
                v_ext = self._data[NN, tstep+1]  # assign extracellular potential value of NN at tstep
            else:
                v_ext = np.zeros(np.shape(self._L[0][gid]))
                for i in range(self._nb_files):
                    period = self._waveforms[i].definition["time"].iloc[-1]
                    simulation_time = (tstep + 1) * sim.dt
                    simulation_time = simulation_time % period
                    v_ext += self._L[i][gid]*self._waveforms[i].calculate(simulation_time)*self._amplitudes[i]

            cell.set_e_extracellular(h.Vector(v_ext))

    def load_comsol(self, comsol_file):
        """Extracts data and headers from comsol.txt. Returns pandas DataFrame.
        The first three columns are the x-, y-, and z-coordinates of the solution nodes.

        :param comsol_file: /path/to/comsol.txt
        :type comsol_file: str
        :return: data
        :rtype: pandas DataFrame
        """
        # Extract column headers and data
        headers = pd.read_csv(comsol_file, sep="\s{3,}", header=None, skiprows=8, nrows=1, engine='python').to_numpy()[0]
        data = pd.read_csv(comsol_file, sep="\s+", header=None, skiprows=9)
        
        # Extract useful info from headers
        headers[0] = headers[0][2:]             # Remove '% ' before first column name
        if headers[3][3] == 'V':                # Convert V to mV if necessary
            data.iloc[:,3:] *= 1000          
        for i,col in enumerate(headers[3:]):    
            if len(data.columns) > 4:   
                headers[i+3] = 1000*float(col[11:])    # Remove superfluous characters 
            else:
                headers[i+3] = 0                  
        
        # Rename column headers
        data.columns = [headers]                
        
        return data

""" 
PSUEDOCODE
----------
__init__:

If comsol_file and waveform are not list
    Make list

If amplitude is not list
    Make list of len(comsol_file)


initialize:

Assert == lens

For i in len(list):
    Load comsol_file
    Create interpolation object
    

step:

For i in len(list):
    ...


"""
