from scipy.interpolate import LinearNDInterpolator as Lip
import sys;
import pandas as pd
import numpy as np
import random

def load_comsol(comsol_file):
        """Extracts data and headers from comsol.txt. Returns pandas DataFrame.
        The first three columns are the x-, y-, and z-coordinates of the solution nodes.
        For a stationary comsol study, the potentials are stored in the fourth column.
        For a time-dependent study, each subsequent column stores the potentials at one timepoints.

        :param comsol_file: (str) "/path/to/comsol.txt"
        :return: (pd DataFrame) Potentials extracted from comsol.txt
        """

        # Extract column headers and data from comsol_file
        headers = pd.read_csv(comsol_file, sep="\s{3,}", header=None, skiprows=8, nrows=1, engine='python')
        headers = headers.to_numpy()[0]
        data = pd.read_csv(comsol_file, sep="\s+", header=None, skiprows=9)

        # Convert V to mV if necessary
        if headers[3][3] == 'V':                        
            data.iloc[:,3:] *= 1000                     

        # Extract useful info from headers
        headers[0] = headers[0][2:]                     # Remove '% ' before first column name
        for i,col in enumerate(headers[3:]):            # Iterate over all elements in the header except first 3
            if len(data.columns) > 4:                   # If time-dependent comsol study
                for j, c in enumerate(col):
                    if c.isdigit():
                        break
                headers[i+3] = 1000*float(col[j:])      # Remove superfluous characters and convert from s to ms
            else:                                       # Else stationary study
                headers[i+3] = 0                        # Rename 4th column
        
        # Rename data with correct column headers
        data.columns = headers

        return data


mask='/Users/ankev/Documents/Github/neural-simulation/v1_Anke/components/stimulation/comsol/mask_small.txt';
data=load_comsol(mask)
points = data[['x', 'y', 'z']]
values = data[0]

interp=Lip(points,values);

def apply_mask_small(positions): # still need to adjust coordinates!! ##
    for coordinate in range(positions.shape[0]):
        interpolated_value=interp(positions[coordinate]);
        #print(interpolated_value)
        if interpolated_value==1000.000000:
            #print(positions[coordinate])
            if positions[coordinate,1]>81:
                positions[coordinate,1]+=1+6.4
            elif positions[coordinate,1]==81:
                positions[coordinate,1]+=random.choice([16.5,-16.5])
            elif positions[coordinate,1]>66 and positions[coordinate,1]<75:
                 positions[coordinate,1]+=1+6.4   
            elif positions[coordinate,1]==66:
                positions[coordinate,1]+=random.choice([7.5,-7.5])
            else:
                positions[coordinate,1]+=-1
            #print(positions[coordinate])
    return positions
        

#positions=np.array([[],[-410,321,870],[6,696,1698],[59.2,499.9,318.9],[-22,509,173],[-22,510,173],[-22,512,173]]);
#mask(positions)
            
    

