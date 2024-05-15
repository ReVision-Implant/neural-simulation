import sys;
module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path);

from file_helper import create_configs;

create_configs('./template.json',4, [0],[10,20], [0], overwrite=True);

