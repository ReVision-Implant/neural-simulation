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
    v1_node_id_slice=v1_node_id[slice_condition];

    #get al information from the original edge file
    target_node_id=hdf5_edges_test.get_target_node_id();
    source_node_id=hdf5_edges_test.get_source_node_id(); #specifies the sender node id of the connection
    edge_group_id=hdf5_edges_test.get_edge_group_id(); # assigns each edge to a specific group of edges
    edge_group_index=hdf5_edges_test.get_edge_group_index();
    edge_type_id=hdf5_edges_test.get_edge_type_id(); #integer to associate an edge to an edge type
    nsyns=hdf5_edges_test.get_nsyns(); #number of connections between the populations
    edges_list=[target_node_id, source_node_id,edge_group_id,edge_group_index,edge_type_id,nsyns];

    #filter the node_ids
    if v1: 
        indices_to_keep=np.isin(target_node_id, v1_node_id_slice) & np.isin(source_node_id,v1_node_id_slice);
    else:
        indices_to_keep=np.isin(target_node_id,v1_node_id_slice);

    filtered_edges_list=[];
    for vector in edges_list:
        filtered_vector=vector[indices_to_keep];
        filtered_edges_list.append(filtered_vector);
    target_node_idf, source_node_idf,edge_group_idf,edge_group_indexf,edge_type_idf,nsynsf=filtered_edges_list;   

    # node id to range is indexed by node id -> not added in new file for now (not necessary?)
    #node_id_to_range_source=hdf5_edges_test.get_node_id_to_range_source();
    #node_id_to_range_sourcef=node_id_to_range_source([slice_condition]);
    #node_id_to_range_target=hdf5_edges_test.get_node_id_to_range_target();
    #node_id_to_range_targetf=node_id_to_range_target([slice_condition]);

    #range_to_edge_id not added for now (not necesarry?)

    #create a new edge file
    import h5py
    import os
    with h5py.File(path_edge_slice, 'w') as slice_edge_test:
        population= os.path.split(path_h5_edgefile)[1][:-12];
        edges=slice_edge_test.create_group('edges');
        v1_to_v1=edges.create_group(population+'_to_v1');
        zero_edge=v1_to_v1.create_group('0');

        v1_to_v1.create_dataset('source_node_id',data=source_node_idf);
        v1_to_v1.create_dataset('target_node_id',data=target_node_idf);
        v1_to_v1.create_dataset('edge_group_id',data=edge_group_idf);
        v1_to_v1.create_dataset('edge_group_index',data=edge_group_indexf);
        v1_to_v1.create_dataset('edge_type_id',data=edge_type_idf);

        zero_edge.create_dataset('nsyns',data=nsynsf);

#test functie
#make_slice_edges('C:/Users/ankev/Documents/GitHub/neural-simulation/v1_Anke/networks_100/network2/v1_nodes.h5','C:/Users/ankev/Documents/GitHub/neural-simulation/v1_Anke/networks_100/network2/lgn_v1_edges.h5','C:/Users/ankev/Documents/GitHub/neural-simulation/v1_Anke/lgn_v1_edges.h5');