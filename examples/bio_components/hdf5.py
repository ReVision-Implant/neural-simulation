import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

class HDF5:

    def __init__(self, dir, plot=False):

        self.dir = dir
        self.file = h5py.File(self.dir, 'r')
        # self.get_positions_v1()
        # self.get_rotations()
        if plot:
            self.plot_positions(labels=['X','Z','Y'])

    def get_positions(self):

        self.name = os.path.split(self.dir)[1][:-9]
        self.positions = self.file['nodes'][self.name]['0']['positions'][:,:]

        self.x_pos = self.positions[:,0]
        self.z_pos = self.positions[:,1]
        self.y_pos = self.positions[:,2]

        return self.positions[:,:]

    def get_positions_v1(self):

        self.name = os.path.split(self.dir)[1][:-9]

        self.x_pos = np.array(self.file['nodes'][self.name]['0']['x'][:])
        self.z_pos = np.array(self.file['nodes'][self.name]['0']['z'][:])
        self.y_pos = np.array(self.file['nodes'][self.name]['0']['y'][:])

        self.positions = np.vstack((self.x_pos, self.y_pos, self.z_pos)).T

        return self.positions

    def get_rotations(self):

        self.rot_x = None
        self.rot_y = None
        self.rot_z = None

        if 'rotation_angle_xaxis' in self.file['nodes'][self.name]['0'].keys():
            self.rot_x = self.file['nodes'][self.name]['0']['rotation_angle_xaxis'][:]
        if 'rotation_angle_yaxis' in self.file['nodes'][self.name]['0'].keys():
            self.rot_y = self.file['nodes'][self.name]['0']['rotation_angle_yaxis'][:]
        if 'rotation_angle_zaxis' in self.file['nodes'][self.name]['0'].keys():
            self.rot_z = self.file['nodes'][self.name]['0']['rotation_angle_zaxis'][:]

        return self.rot_x, self.rot_y, self.rot_z

    def plot_positions(self):
            
            fig, ax = plt.subplots(1,2,figsize=(9,12), subplot_kw={'projection':'3d'})
            
            self.subplot(ax[0], 'YZ')
            ax[0].view_init(elev=0, azim=0)

            self.subplot(ax[1], 'XY')
            ax[1].view_init(elev=90., azim=270)
            plt.show()

    def subplot(self, ax, title):

        ax.scatter(self.x_pos,self.y_pos,self.z_pos, marker='.', s=3)
        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')


# nodes = HDF5('/users/students/r0754386/Documents/bmtk/examples/v1/networks_rebuilt/network/v1_nodes.h5')
# print(nodes.get_positions_v1())