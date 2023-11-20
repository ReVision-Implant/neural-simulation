def make_slice_edges(path_h5file,path_h5_edgefile,path_edge_slice,v1=False):
    #remark: for the edges of the slice it is important to discriminate between connection within v1 vs with v1
    import matplotlib.pyplot as plt;
    import numpy as np;
    from hdf5_edge import HDF5_EDGE;
    hdf5_edges_test=HDF5_EDGE(path_h5_edgefile);

    import sys;
    module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox/hdf5.py';
    sys.path.append(module_path);
    from hdf5 import HDF5;
    hdf5_test=HDF5(path_h5file,v1=True);
    positions=hdf5_test.get_positions(v1=True);
    slice_condition=(positions[:,0] >= -150) & (positions[:,0]<=150);

    #get node_id from node file
    v1_node_id=hdf5_test.get_node_id();
    #print(v1_node_id[:20]);
    v1_node_id_slice=v1_node_id[slice_condition];
    #print(v1_node_id_slice[:20]);
    new_id=np.arange(len(v1_node_id_slice));
    matrix_node_id=np.column_stack((v1_node_id_slice, new_id));
    

    #get al information from the original edge file
    target_node_id=hdf5_edges_test.get_target_node_id();
    #print(target_node_id[3000:3063]);
    source_node_id=hdf5_edges_test.get_source_node_id(); #specifies the sender node id of the connection
    edge_group_id=hdf5_edges_test.get_edge_group_id(); # assigns each edge to a specific group of edges
    edge_group_index=hdf5_edges_test.get_edge_group_index();
    edge_type_id=hdf5_edges_test.get_edge_type_id(); #integer to associate an edge to an edge type
    nsyns=hdf5_edges_test.get_nsyns(); #number of connections between the populations
   

    #filter the node_ids
    if v1: 
        indices_to_keep=np.isin(target_node_id, v1_node_id_slice) & np.isin(source_node_id,v1_node_id_slice);    
    else:
        indices_to_keep=np.isin(target_node_id,v1_node_id_slice);    

    #print(indices_to_keep[:20])
    target_node_idf=target_node_id[indices_to_keep];
    source_node_idf=source_node_id[indices_to_keep];
    edge_group_idf=edge_group_id[indices_to_keep];
    edge_group_indexf=edge_group_index[indices_to_keep];
    edge_type_idf=edge_type_id[indices_to_keep];
    nsynsf=nsyns[indices_to_keep];

    #some debugging code
    #print(target_node_idf[:20])
    #index=np.where(matrix_node_id[:,0]==target_node_idf[0])[0];
    #index_value=index[0];
    #print(index_value)
    
    prev_value= None
    prev_index= None

    for node_id, value_to_find in enumerate (target_node_idf):
        if value_to_find != prev_value:
            index=np.where(matrix_node_id[:,0]==value_to_find)[0]; #returns an array with one value
            index_value=index[0]; #gets the value of that array
            prev_value=value_to_find
            prev_index=index_value

        else:
            index_value=prev_index

        target_node_idf[node_id]=matrix_node_id[index_value,1];

    print('the shape of the target node id  is',target_node_idf.shape,'the first ten values are',target_node_idf[:10]);   

    if v1:
        prev_value_s= None
        prev_index_s= None
        for node_id, value_to_find in enumerate (source_node_idf):
            if value_to_find != prev_value_s:
                index=np.where(matrix_node_id[:,0]==value_to_find)[0]; #returns an array with one value
                index_value=index[0]; #gets the value of that array
                prev_value_s=value_to_find
                prev_index_s=index_value
            else:
                index_value=prev_index_s
                
            source_node_idf[node_id]=matrix_node_id[index_value,1];
        print('the shape of the source node id is',source_node_idf.shape,'the first ten values are',source_node_idf[:10]);

    #create a new edge file
    import h5py
    import os
    with h5py.File(path_edge_slice, 'w') as slice_edge_test:
        population= os.path.split(path_h5_edgefile)[1][:-12];
        edges=slice_edge_test.create_group('edges');
        population_to_v1=edges.create_group(population+'_to_v1');
        zero_edge=population_to_v1.create_group('0');

        population_to_v1.create_dataset('source_node_id',data=source_node_idf);
        population_to_v1.create_dataset('target_node_id',data=target_node_idf);
        population_to_v1.create_dataset('edge_group_id',data=edge_group_idf);
        population_to_v1.create_dataset('edge_group_index',data=edge_group_indexf);
        population_to_v1.create_dataset('edge_type_id',data=edge_type_idf);

        zero_edge.create_dataset('nsyns',data=nsynsf);

        #add attributes to the source_node_id and target_node_id
        source_node_id_dataset=population_to_v1['source_node_id'];
        source_node_id_dataset.attrs['node_population']=population;
        
        target_node_dataset=population_to_v1['target_node_id'];
        target_node_dataset.attrs['node_population']='v1';

#test functie
#make_slice_edges('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_1/v1_nodes.h5','/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_1/v1_v1_edges.h5','/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/mouse_1/v1_v1_edges.h5',v1=True);