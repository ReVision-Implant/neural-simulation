import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

class HDR5:

    def __init__(self, dir):

        self.dir = dir
        self.file = h5py.File(self.dir, 'r')
        
    def get_positions(self, plot=False):

        self.name = os.path.split(self.dir)[1][:-9]
        self.positions = self.file['nodes'][self.name]['0']['positions'][:,:]
        self.x_pos = self.positions[:,0]
        self.z_pos = self.positions[:,1]
        self.y_pos = self.positions[:,2]

        if plot:

            fig, ax = plt.subplots(1,2,figsize=(9,12), subplot_kw={'projection':'3d'})
            self.plot_positions(ax[0], labels=['X','Z','Y'])
            ax[0].set_title('Side view')
            ax[0].view_init(elev=5., azim=0)

            self.plot_positions(ax[1], labels=['X','Z','Y'])
            ax[1].set_title('Bird\'s eye view')
            ax[1].view_init(elev=90., azim=0)
            plt.show()

        return self.positions[:,:]

    def plot_positions(self, ax, labels):
        ax.scatter(self.x_pos,self.y_pos,self.z_pos, marker='.', s=3)

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])

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
