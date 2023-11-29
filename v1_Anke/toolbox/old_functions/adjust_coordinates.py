import h5py
import numpy as np
import shutil
import sys;
module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path);
from hdf5 import HDF5

def adjust_coordinates(coordinates, volumes, displacement=0.1e-6):
    adjusted_coordinates = np.copy(coordinates)

    for volume in volumes:
        min_values = volume['min']
        max_values = volume['max']

        # Extract middle y-coordinate of the volume
        middle_y = (min_values[1] + max_values[1]) / 2.0

        # Check if each coordinate falls within the volume
        within_volume_mask = np.all((min_values <= adjusted_coordinates) & (adjusted_coordinates <= max_values), axis=1)

        # Adjust y-coordinate based on its relationship to the middle y-coordinate
        left_of_middle_mask = adjusted_coordinates[within_volume_mask, 1] < middle_y
        right_of_middle_mask = adjusted_coordinates[within_volume_mask, 1] > middle_y
        in_the_middle_mask = ~left_of_middle_mask & ~right_of_middle_mask

        # Randomly decide whether to move left or right for coordinates in the middle
        if np.any(in_the_middle_mask):
            move_left_mask = in_the_middle_mask & (np.random.rand(np.sum(in_the_middle_mask)) < 0.5)
        else:
            move_left_mask = np.array([], dtype=bool)

        # Ensure all masks are defined
        if not np.any(left_of_middle_mask):
            left_of_middle_mask = np.array([], dtype=bool)
        if not np.any(right_of_middle_mask):
            right_of_middle_mask = np.array([], dtype=bool)
        if not np.any(in_the_middle_mask):
            in_the_middle_mask = np.array([], dtype=bool)

        # Adjust coordinates
        adjusted_coordinates[within_volume_mask, 1] += np.where(left_of_middle_mask, -displacement, 0)
        adjusted_coordinates[within_volume_mask, 1] += np.where(right_of_middle_mask, displacement, 0)
        adjusted_coordinates[within_volume_mask, 1] += np.where(move_left_mask, -displacement, 0)


    return adjusted_coordinates

def change_coordin(input_file, output_path,volumes):
    # get coordinates original file
    file=HDF5(input_file,v1=True);
    coordinates=file.get_positions(v1=True);

    # Adjust coordinates
    adjusted_coordinates = adjust_coordinates(coordinates, volumes)

    #copy the file to new location
    shutil.copy2(input_file, output_path)
    
    #delete old coordinates and insert new ones
    with h5py.File(output_path,'a') as file:
        zero=file['nodes/v1/0']
        del zero['x'];
        del zero['y'];
        del zero['z'];

        zero.create_dataset('x',data=adjusted_coordinates[:,0]);
        zero.create_dataset('y',data=adjusted_coordinates[:,1]);
        zero.create_dataset('z',data=adjusted_coordinates[:,2]);
        