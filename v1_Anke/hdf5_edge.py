import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

class HDF5_EDGE:
    """ Helper class to extract edge data from .h5 files.
    """
    def __init__(self, dir):
        """Read .h5 file and call main functions.

        :param dir: path to .h5 file.
        :type dir: path
        """        
        self.dir = dir
        self.file = h5py.File(self.dir, 'r')
    
    def get_target_node_id(self):

        self.name = os.path.split(self.dir)[1][:-12]
        target_node_id=np.array(self.file['edges'][self.name+'_to_v1']['target_node_id'][:])
        return target_node_id
    
    def get_source_node_id(self):
        self.name = os.path.split(self.dir)[1][:-12]
        source_node_id=np.array(self.file['edges'][self.name+'_to_v1']['source_node_id'][:])
        return source_node_id
    
    def get_edge_group_id(self):
        self.name = os.path.split(self.dir)[1][:-12]
        edge_group_id=np.array(self.file['edges'][self.name+'_to_v1']['edge_group_id'][:])
        return edge_group_id
    
    def get_edge_group_index(self):
        self.name = os.path.split(self.dir)[1][:-12]
        edge_group_index=np.array(self.file['edges'][self.name+'_to_v1']['edge_group_index'][:])
        return edge_group_index
    
    def get_edge_type_id(self):
        self.name = os.path.split(self.dir)[1][:-12]
        edge_type_id=np.array(self.file['edges'][self.name+'_to_v1']['edge_type_id'][:])
        return edge_type_id
    
    def get_nsyns(self):
        self.name = os.path.split(self.dir)[1][:-12]
        nsyns=np.array(self.file['edges'][self.name+'_to_v1']['0']['nsyns'][:])
        return nsyns
    
    #def get_node_id_to_range_source(self):
        node_id_to_range_source=np.array(self.file['edges']['v1_to_v1']['indices']['source_to_target']['node_id_to_range'][:])
        return node_id_to_range_source
    
    #def get_range_to_edge_id_source(self):
        range_to_edge_id_source=np.array(self.file['edges']['v1_to_v1']['indices']['source_to_target']['range_to_edge_id'][:])
        return range_to_edge_id_source
    
    #def get_node_id_to_range_target(self):
        node_id_to_range_target=np.array(self.file['edges']['v1_to_v1']['indices']['target_to_source']['node_id_to_range'][:])
        return node_id_to_range_target
    
    #def get_range_to_edge_id_target(self):
        range_to_edge_id_target=np.array(self.file['edges']['v1_to_v1']['indices']['target_to_source']['range_to_edge_id'][:])
        return range_to_edge_id_target
    