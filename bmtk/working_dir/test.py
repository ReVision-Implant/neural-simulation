""" COMSOL -> NEURON link
IN: 
    - COMSOL: 
        - file containing positions and extracellular potential distribution
        - mesh file
    - NEURON:
        - position of NEURON node (=compartment)

OUTPUT: (interpolated/nearest neighbour) extracellular potential at NEURON node
"""

import pandas as pd
import numpy as np
from scipy.interpolate import NearestNDInterpolator as NN
COMSOL = pd.read_csv("Untitled.txt", sep="\s+", header=8, usecols=[0,1,2,3], names=['x','y','z','V'])

interp = NN(COMSOL[['x','y','z']], COMSOL['V'])
r05 = 500*np.random.rand(10000,3)
v_ext = np.expand_dims(interp(r05),1)
print(np.hstack((r05[:25],v_ext[:25])))

def get_vext(self):
    interp = NN(COMSOL[['x','y','z']], COMSOL['V'])
    r05 = self.seg_coords.p05
    nseg = r05.shape[1]
    v_extracellular = interp(r05)
    vext_vec = h.Vector(v_extracellular)
    return vext_vec
