#helperfunction to make a slice from an h5 file
def make_slice_nodes(path_h5file,path_slice):
    import sys;
    module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox/hdf5.py';
    sys.path.append(module_path);

    import matplotlib.pyplot as plt;
    import numpy as np;
    from hdf5 import HDF5;
    import shutil
    import h5py

    hdf5_test=HDF5(path_h5file+'v1_nodes.h5',v1=True);
    #hdf5_test.plot_positions()

    #get positions that lie within slice
    positions=hdf5_test.get_positions(v1=True);
    slice_condition=(positions[:,0] >= -150) & (positions[:,0]<=150);
    slice_positions=positions[slice_condition];

    #get rotations in the slice
    rot_x,rot_y,rot_z=hdf5_test.get_rotations()
    #print(rot_x.shape,rot_y.shape,rot_z.shape)
    rotations=np.vstack((rot_x,rot_y,rot_z)).T;
    slice_rotations=rotations[slice_condition];

    #get the tuning angle of the positions within the slice
    tuning_angles_original=hdf5_test.get_tuning_angles();
    tuning_angles_slice=tuning_angles_original[slice_condition];

    #get the group_node_ids within the slice
    group_node_ids=hdf5_test.get_node_group_ids();
    group_node_ids_slice=group_node_ids[slice_condition];

    #get the node_group_index
    node_group_index=hdf5_test.get_node_group_index();
    node_group_index_slice=node_group_index[slice_condition];

    #get node_id
    node_id=hdf5_test.get_node_id();
    node_id_slice=node_id[slice_condition];
    #print(node_id_slice.shape)

    #get node_type_id
    node_type_id=hdf5_test.get_node_type_id();
    node_type_id_slice=node_type_id[slice_condition];

    #copy the original HDF5 file to the new location
    shutil.copyfile(path_h5file+'v1_nodes.h5',path_slice+'v1_nodes.h5')

    #open the copied file for modification
    with h5py.File(path_slice+'v1_nodes.h5', 'a') as slice_test:
        if 'nodes' in slice_test: # deletes the information from prev file 
            del slice_test['nodes']
        nodes=slice_test.create_group('nodes');
        v1=nodes.create_group('v1');
        zero=v1.create_group('0');

        zero.create_dataset('x', data=slice_positions[:,0]);
        zero.create_dataset('y', data=slice_positions[:,1]);
        zero.create_dataset('z', data=slice_positions[:,2]);

        #add rotations
        zero.create_dataset('rotation_angle_xaxis', data=slice_rotations[:,0]);
        zero.create_dataset('rotations_angle_yaxis',data=slice_rotations[:,1]);
        zero.create_dataset('rotations_angle_zaxis',data=slice_rotations[:,2]);

        #add tuning angle
        zero.create_dataset('tuning_angle',data=tuning_angles_slice);

        #add group node id: assigns each node to a specific group of nodes
        v1.create_dataset('node_group_id',data=group_node_ids_slice);

        #add group node index: indicates the index within that group that contains all the attributes for a particular node
        v1.create_dataset('node_group_index',data=node_group_index_slice);

        #add node_id: assigns a key to uniquely identify and lookup a node within a population
        v1.create_dataset('node_id',data=node_id_slice);

        #add node_type_id
        v1.create_dataset('node_type_id',data=node_type_id_slice);

#create a new csv file -> does not seem to be necessary 
    #import pandas as pd

    #load the csv file into a pandas dataframe 
    #data=pd.read_csv(path_h5file+'v1_node_types.csv')
    #data=data['node_type_id pop_name ei model_type population model_processing morphology dynamics_params location model_template'].str.split(' ', expand=True);
    #print(data.columns)

    #filter the data frame
    #filtered_data=data[~data['node_type_id'].sin(node_type_id_slice)]

    #save to a new csv file
    #filtered_data.to_csv(path_slice+'v1_node_types.csv', index=False);
 
