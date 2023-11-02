def slice_mouse(mouse_number):
    import sys;
    module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke';
    sys.path.append(module_path);

    from make_slice_nodes import make_slice_nodes;
    from make_slice_edges import make_slice_edges;

    mouse_number=str(mouse_number);

    make_slice_nodes('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/v1_nodes.h5', '/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/slice_mouse_'+mouse_number+'/v1_nodes.h5');
    make_slice_edges('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/v1_nodes.h5', '/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/v1_v1_edges.h5','/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/slice_mouse_'+mouse_number+'/v1_v1_edges.h5',v1=True);
    make_slice_edges('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/v1_nodes.h5', '/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/lgn_v1_edges.h5','/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/slice_mouse_'+mouse_number+'/lgn_v1_edges.h5',v1=False);
    make_slice_edges('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/v1_nodes.h5', '/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/virtual_mice/mouse_'+mouse_number+'/bkg_v1_edges.h5','/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/slices/slice_mouse_'+mouse_number+'/bkg_v1_edges.h5',v1=False);

#copy the lgn, bkg nodes and csv files manually 