import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
import shutil

class HDF5:
    """Helper class to extract data from .h5 files.
    """
    def __init__(self, dir, v1=False, plot=False):
        """Read .h5 file and call main functions.

        :param dir: path to .h5 file.
        :type dir: path
        :param v1: set to True if nodes.h5 represents V1 column, defaults to False.
        :type v1: bool, optional
        :param plot: if True, show plot of node positions. defaults to False.
        :type plot: bool, optional
        """        
        self.dir = dir
        self.file = h5py.File(self.dir, 'r+')
        #self.v1 = v1
        self.get_positions(v1)
        rot = self.get_rotations()
        
        if plot:
            self.plot_positions(labels=['X','Z','Y'])
        
    def get_positions(self,v1=False):
        """Get node positions.

        :param v1: if True, interpret nodes.h5 as V1 column. defaults to False.
        :type v1: bool, optional
        :return: node positions
        :rtype: ndarray
        """
        self.name = os.path.split(self.dir)[1][:-9]
        
        if v1:

            self.x_pos = np.array(self.file['nodes'][self.name]['0']['x'][:])
            self.y_pos = np.array(self.file['nodes'][self.name]['0']['y'][:])
            self.z_pos = np.array(self.file['nodes'][self.name]['0']['z'][:])
            self.positions = np.vstack((self.x_pos, self.y_pos, self.z_pos)).T

            return self.positions
            
        else:

            self.positions = self.file['nodes'][self.name]['0']['positions'][:,:]
            self.x_pos = self.positions[:,0]
            self.y_pos = self.positions[:,1]
            self.z_pos = self.positions[:,2]

            return self.positions[:,:]

    def get_rotations(self):
        """Get node rotations.

        :return: node rotations
        :rtype: ndarray
        """
        ### Initialise
        self.x_rot = None
        self.y_rot = None
        self.z_rot = None

        ### If rotation angles exist, set as attributes
        if 'rotation_angle_xaxis' in self.file['nodes'][self.name]['0'].keys():
            self.x_rot = self.file['nodes'][self.name]['0']['rotation_angle_xaxis'][:]
        if 'rotation_angle_yaxis' in self.file['nodes'][self.name]['0'].keys():
            self.y_rot = self.file['nodes'][self.name]['0']['rotation_angle_yaxis'][:]
        if 'rotation_angle_zaxis' in self.file['nodes'][self.name]['0'].keys():
            self.z_rot = self.file['nodes'][self.name]['0']['rotation_angle_zaxis'][:]
        
        ### Create dictionary of present rotation angles {'axis_rot':self.axis_rot}
        self.rotations = {a:getattr(self, a) for a in ['x_rot', 'y_rot', 'z_rot'] if hasattr(self,a)}
        self.rotations = {k:v for k, v in self.rotations.items() if v is not None}

        return self.x_rot,self.y_rot,self.z_rot

    def plot_positions(self):
        """Plot node positions.
        """
        ### Create 1x2 subplot    
        _, ax = plt.subplots(1,2,figsize=(9,12), subplot_kw={'projection':'3d'})
        
        ### Plot first perspective
        self.subplot(ax[0], 'YZ')
        ax[0].view_init(elev=0, azim=0)

        ### Plot second perspective
        self.subplot(ax[1], 'XY')
        ax[1].view_init(elev=90., azim=270)
        plt.show()

    def subplot(self, ax, title):
        """Helper function for plot_position()

        :param ax: matplotlib axes object.
        :type ax: matplotlib.axes object
        :param title: subplot title.
        :type title: str
        """
        ### Set attributes
        ax.scatter(self.x_pos,self.y_pos,self.z_pos, marker='.', s=3)
        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    def convert_pos_to_v1(self):
        """Convert 'positions' dataset into three datasets 'x', 'y', 'z' in nodes/name/0 for use in VND.
        """
        ### If nodes are already in v1 format, do nothing
        if self.v1:
            pass
        ### Else, try creating new datasets 'x', 'y', 'z' or overwrite if they already exist
        else: 
            try:
                self.file.create_dataset('nodes/'+self.name+'/0/x', data=self.x_pos)
                self.file.create_dataset('nodes/'+self.name+'/0/y', data=self.y_pos)
                self.file.create_dataset('nodes/'+self.name+'/0/z', data=self.z_pos)
            except:
                self.file['nodes'][self.name]['0']['x'][...] = self.x_pos
                self.file['nodes'][self.name]['0']['y'][...] = self.y_pos
                self.file['nodes'][self.name]['0']['z'][...] = self.z_pos
            self.file.close()
    
        return
    
    #function for getting the tuning angles
    def get_tuning_angles(self):
        self.name = os.path.split(self.dir)[1][:-9]
        tuning_angles=np.array(self.file['nodes'][self.name]['0']['tuning_angle'][:])
        return tuning_angles
    
    #function to get node_group_ids
    def get_node_group_ids(self):
        self.name=os.path.split(self.dir)[1][:-9]
        node_group_id=np.array(self.file['nodes'][self.name]['node_group_id'][:])
        return node_group_id

    #function to get node_group_index
    def get_node_group_index(self):
        self.name=os.path.split(self.dir)[1][:-9]
        node_group_index=np.array(self.file['nodes'][self.name]['node_group_index'][:])
        return node_group_index   

    #function to get node_id
    def get_node_id(self):
        self.name=os.path.split(self.dir)[1][:-9]
        node_id=np.array(self.file['nodes'][self.name]['node_id'][:])
        return node_id
    
    #function to get node_type_id
    def get_node_type_id(self):
        self.name=os.path.split(self.dir)[1][:-9]
        node_type_id=np.array(self.file['nodes'][self.name]['node_type_id'][:])
        return node_type_id

class SpikeScaler():
    """Scale spike timings for use in VND.
    """
    def __init__(self, path, scale):
        """Read .h5 file and call main function.

        :param path: Path to .h5 file.
        :type path: path
        :param scale: Scale factor for the spike timings. The time axis is multiplied by this scale factor.
        :type scale: int
        """
        ### Set attributes
        self.source = os.path.abspath(__file__) + '/../../../' + path
        self.dest = os.path.dirname(self.source)+"/spikes_scaled.h5"
        self.scale = scale
        self.scale_spike_timings()

    def scale_spike_timings(self):
        """Copies the original spikes.h5 and changes the time scale in the copy.
        """
        ### Copy original file into scaled file
        self.orig_file = h5py.File(self.source, 'r')
        self.scaled_file = h5py.File(self.dest, 'w')
        self.scaled_file.copy(self.orig_file['spikes'], 'spikes')
        
        ### Scale timestamps 
        for i in list(self.scaled_file['spikes'].keys()):
            timestamps = self.scaled_file['spikes/'+i+'/timestamps'][:]
            del self.scaled_file['spikes/'+i+'/timestamps']
            try:
                self.scaled_file.create_dataset('spikes/'+i+'/timestamps', data=timestamps*self.scale)
            except:
                self.scaled_file['spikes/'+i+'/timestamps'][:] = timestamps*self.scale
            units = self.orig_file['spikes/'+i+'/timestamps'].attrs.__getitem__('units')
            self.scaled_file['spikes/'+i+'/timestamps'].attrs.__setitem__('units', units)

        ### Close both files
        self.orig_file.close()
        self.scaled_file.close()
        
        return

#if __name__ == "__main__":
    #root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    #hf = HDF5(root + '/examples/comsol/network/column_nodes.h5')
    #hf.convert_pos_to_v1()

    #spikes = SpikeScaler('/examples/comsol/outputs/output_stat/spikes.h5', 100)
