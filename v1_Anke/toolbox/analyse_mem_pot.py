from bmtk.analyzer.compartment import plot_traces
import sys;
module_path='/scratch/leuven/356/vsc35693/neural-simulation/';
sys.path.append(module_path)

dir='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_1/config/bkg/bkg_membrane_pot/exp_C/config.exp1_bkg_vm_C'
title_plot='background membrane potential v1 model'

plot_traces(config_file=dir+'.json', report_name='v_report', title=title_plot, save_as= '/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_1/output/bkg/bkg_membrane_pot/exp_C/plot_vm.png')
