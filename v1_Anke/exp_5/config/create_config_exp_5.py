import sys;
module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path);

from file_helper import create_configs;

create_configs('./template.json',5, [0,1,2,3,4,5,6,7,8,9],[10,20],[2], overwrite=True);

