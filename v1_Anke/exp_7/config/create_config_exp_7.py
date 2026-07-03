import sys;
#module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
module_path=module_path = r"C:\Users\ankev\OneDrive\Documenten\Github\ReVision\neural-simulation\v1_Anke\toolbox"
sys.path.append(module_path);

from file_helper import create_configs;

create_configs('./template.json',7, [0,1,2,3],[10,20],[0,1,2], overwrite=True);

