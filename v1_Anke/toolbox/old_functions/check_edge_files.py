def check_edge_files(path_nodes_file,path_edge_file):
    import sys;
    module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
    sys.path.append(module_path);

    from hdf5 import HDF5;
    hdf5_test=HDF5(path_nodes_file,v1=True);

    from hdf5_edge import HDF5_EDGE;
    hdf5_edges_test=HDF5_EDGE(path_edge_file)

    v1_node_id=hdf5_test.get_node_id()
    source_node_id=hdf5_edges_test.get_source_node_id()

    set_source_node_id=set(source_node_id);
    set_v1_node_id=set(v1_node_id);

    if set_source_node_id.issubset(set_v1_node_id):
        print('ok')
    else:
        missing_nodes=set_source_node_id.difference(set_v1_node_id)
        print('following nodes are missing in v1_node_id:',missing_nodes);

#check_edge_files('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/mouse_0/v1_nodes.h5','/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/mouse_0/v1_v1_edges.h5')

#check the population names
import sonata
edge_file=sonata.file(data_files='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/mouse_0/v1_nodes.h5');
population_names=edge_file.populations.keys()

for name in population_names:
    print(name)