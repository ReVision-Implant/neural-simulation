import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

class HDF5:
    """ Helper class to extract spike data from output.h5 files.
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
        self.file = h5py.File(self.dir, 'r')
        self.get_positions(v1)
        self.get_rotations()
        if plot:
            self.plot_positions(labels=['X','Z','Y'])

    def get_positions(self, v1=False):
        """Get node positions.

        :param v1: if True, interpret nodes.h5 as V1 column. defaults to False.
        :type v1: bool, optional
        :return: node positions
        :rtype: ndarray
        """
        self.name = os.path.split(self.dir)[1][:-9]
        
        if v1:

            self.x_pos = np.array(self.file['nodes'][self.name]['0']['x'][:])
            self.z_pos = np.array(self.file['nodes'][self.name]['0']['z'][:])
            self.y_pos = np.array(self.file['nodes'][self.name]['0']['y'][:])
            self.positions = np.vstack((self.x_pos, self.y_pos, self.z_pos)).T

            return self.positions
            
        else:

            self.positions = self.file['nodes'][self.name]['0']['positions'][:,:]
            self.x_pos = self.positions[:,0]
            self.z_pos = self.positions[:,1]
            self.y_pos = self.positions[:,2]

            return self.positions[:,:]

    def get_rotations(self):
        """Get node rotations.

        :return: node rotations
        :rtype: ndarray
        """
        self.x_rot = None
        self.y_rot = None
        self.z_rot = None

        if 'rotation_angle_xaxis' in self.file['nodes'][self.name]['0'].keys():
            self.x_rot = self.file['nodes'][self.name]['0']['rotation_angle_xaxis'][:]
        if 'rotation_angle_yaxis' in self.file['nodes'][self.name]['0'].keys():
            self.y_rot = self.file['nodes'][self.name]['0']['rotation_angle_yaxis'][:]
        if 'rotation_angle_zaxis' in self.file['nodes'][self.name]['0'].keys():
            self.z_rot = self.file['nodes'][self.name]['0']['rotation_angle_zaxis'][:]
        
        self.rotations = np.vstack((self.x_rot, self.y_rot, self.z_rot)).T

        return self.rotations

    def plot_positions(self):
        """Plot node positions.
        """            
        _, ax = plt.subplots(1,2,figsize=(9,12), subplot_kw={'projection':'3d'})
        
        self.subplot(ax[0], 'YZ')
        ax[0].view_init(elev=0, azim=0)

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
        ax.scatter(self.x_pos,self.y_pos,self.z_pos, marker='.', s=3)
        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')