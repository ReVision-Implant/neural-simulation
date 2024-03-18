import json
import os
import pandas as pd

def create_configs(template, exp, patterns, amplitudes, mice, overwrite=False):
    """Create config files according to specified stimulation parameters. 
    All input parameters will automatically be converted to <parameter>_<name>.

    :param template: Template config.json that is modified to generate all other files. 
        For each experiment, it is generally a good idea to manually configure one config.json file.
        If that works, this function will then do the rest. 
    :param exp: Experiment name. int/str.
    :param patterns: Pattern names. int/str or list thereof.
    :param amplitudes: Amplitudes in uA. int/str or list thereof.
    :param mice: Mouse names. int/str or list thereof.
    :param overwrite: Overwrite existing files if True. defaults to False.
    :return: None
    
    | Example:
    | Passing arguments 1, '1', [1], ['1'], 'pattern1', ['pattern1'], 'pattern_1', ['pattern_1'] will all be formatted to ['pattern_1']
    """

    # Convert all parameters into [name]_[value] format
    exp, patterns, amplitudes, mice = format_params(exp, patterns, amplitudes, mice)

    # Load template config.json
    with open(template, 'r') as read_file:
        data = json.load(read_file)

    # Nested loop over all specified parameters (converted to str)
    for pattern in patterns:
        
        # Get pattern parameters from pattern.csv
        stim_params = get_stim_params(exp, pattern)

        # Set paths to comsol.txt and waveform.csv files in config.json
        comsol_files = [os.path.join('$STIM_DIR', 'comsol', file) for file in stim_params['comsol_files']]
        waveforms = [os.path.join('$STIM_DIR', 'waveforms', file) for file in stim_params['waveforms']]
        data['inputs']['Extracellular_Stim']['comsol_files'] = comsol_files
        data['inputs']['Extracellular_Stim']['waveforms'] = waveforms

        for amplitude in amplitudes:
            
            # Set amplitude in config.json
            amplitude_values = [i*float(amplitude[10:]) for i in stim_params['amplitudes']]
            data['inputs']['Extracellular_Stim']['amplitudes'] = amplitude_values

            for mouse in mice:
                
                # Set mouse dir in config.json
                network_dir = os.path.join(os.path.dirname(data['manifest']['$NETWORK_DIR']), mouse)
                data['manifest']['$NETWORK_DIR'] = network_dir
                                            
                # Set output dir in config.json
                output_dir = os.path.join('$BASE_DIR', exp , 'output', pattern, amplitude, mouse)
                data['manifest']["$OUTPUT_DIR"] = output_dir

                # Get path for config.json and create directory if it does not exist
                root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                file_dir = os.path.join(root, exp, 'config', pattern, amplitude)
                if not os.path.exists(file_dir):
                    os.makedirs(file_dir)
                
                # Save config.json unless it already exists and overwrite is False
                file_name =  'config.' + '-'.join((exp[0]+exp[3:], pattern[0]+pattern[7:], amplitude[0]+amplitude[9:], mouse[0]+mouse[5:])) + '.json'
                if (overwrite == True) or (not os.path.exists(os.path.join(file_dir, file_name))):
                    with open(os.path.join(file_dir, file_name), 'w') as write_file:
                        json.dump(data, write_file, indent=2)
                        print('created file: ',file_name, 'in', file_dir)
    return

def delete_configs(exp, patterns=None, amplitudes=None, mice=None):
    """Delete config files according to specified stimulation parameters. Also deletes empty folders.
    All input parameters will automatically be converted to <parameter>_<name>.
    
    :param exp: Experiment name. int/str.
    :param patterns: Pattern names. int/str or list thereof.
    :param amplitudes: Amplitudes in uA. int/str or list thereof.
    :param mice: Mouse names. int/str or list thereof.
    :return: None
    
    Example::
    
        Passing arguments 1, '1', [1], ['1'], 'pattern1', ['pattern1'], 'pattern_1', ['pattern_1'] will all be formatted to ['pattern_1']
    """

    # Convert all parameters into [name]_[value] format
    exp, patterns, amplitudes, mice = format_params(exp, patterns, amplitudes, mice)

    # Get directory where all exp_? folders are
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 

    # If any of the input arguments is None, skip this statement and just remove empty directories
    if patterns is None or amplitudes is None or mice is None:
        pass
    else:
        # Nested loop over all specified parameters (converted to strings)
        for pattern in patterns:
            for amplitude in amplitudes:
                for mouse in mice:
                    
                    # Get file name and path 
                    file_dir = os.path.join(root, exp, 'config', pattern, amplitude)
                    file_name =  'config.' + '-'.join((exp[0]+exp[3:], pattern[0]+pattern[7:], amplitude[0]+amplitude[9:], mouse[0]+mouse[5:])) + '.json'
                    
                    # Remove file if it exists
                    if os.path.exists(os.path.join(file_dir, file_name)):
                        os.remove(os.path.join(file_dir, file_name))
                        print('removed file: ', file_name)

    # Remove empty directories
    walk = list(os.walk(os.path.join(root,exp,'config')))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)
            print('removed folder: ', path)

    # Delete exp folder is empty if empty
    path = os.path.join(root,exp)
    if len(os.listdir(path)) == 0:
                os.rmdir(path)
                print('removed folder: ', path)

    return

def get_stim_params(exp, pattern):
    """Gets extracellular stimulation parameters from the specified pattern.csv in ./stimulation/patterns/exp/.
    These are used to write in the inputs/Extracellular_Stim section of config.json. 

    :param exp: Experiment name. 
    :param pattern: Pattern name. 
    :return: Dictionary containing all stimulation parameters
    """    

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pattern = str(pattern)+'.csv' if str(pattern)[-4:] != '.csv' else pattern
    path = os.path.join(root,'components', 'stimulation', 'patterns', exp, pattern)
    df = pd.read_csv(path, sep='\s+')
    
    comsol_files = df['comsol'].tolist()
    waveforms = df['waveform'].tolist()
    amplitudes = df['amplitude'].tolist()
    
    return dict(comsol_files=comsol_files, waveforms=waveforms, amplitudes=amplitudes)

def format_params(exp, patterns, amplitudes, mice):
    """Standardise parameters to <parameter>_<name>.
    Example:: 1, '1', [1], ['1'], 'pattern1', ['pattern1'], 'pattern_1', ['pattern_1'] will all be formatted to ['pattern_1']

    :param exp: Experiment name. int/str.
    :param patterns: Pattern names. int/str or list thereof.
    :param amplitudes: Amplitudes in uA. int/str or list thereof.
    :param mice: Mouse names. int/str or list thereof.
    :return: standardised parameter names
    """    
    exp = [exp] if not isinstance(exp, list) else exp
    exp = [format_param(str(i), 'exp') for i in exp]
    exp = exp[0] if len(exp) == 1 else exp

    patterns = [patterns] if not isinstance(patterns, list) else patterns
    patterns = [format_param(str(i), 'pattern') for i in patterns]

    amplitudes = [amplitudes] if not isinstance(amplitudes, list) else amplitudes
    amplitudes = [format_param(str(i), 'amplitude') for i in amplitudes]

    mice = [mice] if not isinstance(mice, list) else mice
    mice = [format_param(str(i), 'mouse') for i in mice]

    return exp, patterns, amplitudes, mice

def format_param(input, name):
    """Format parameter to <parameter>_<name>.

    Example::

        '1', 'pattern1', 'pattern_1', will all be formatted to 'pattern_1'

    :param input: _description_
    :param name: _description_
    :return: _description_
    """    
    n = len(name)
    if input[:n+1] == name + '_':
        output = input
    elif input[:n] == name:
        output = name + '_' + input[n:]
    else:
        output = name + '_' + input

    return output

def get_dirs(exp, patterns, amplitudes, mice):
    """Gets the directories to the nodes.h5, spikes.csv, spikes_bkg.csv, pattern.csv, and electrodes.csv files.
    This function does this for all combinations of the input arguments and returns them in lists.

    :param exp: Experiment name. int/str.
    :param patterns: Pattern names. int/str or list thereof.
    :param amplitudes: Amplitudes in uA. int/str or list thereof.
    :param mice: Mouse names. int/str or list thereof.
    :return: Dictionary containing where the values are lists containing the dirs corresponding the key.
    """    
    # Initialisation
    nodes_dirs = []
    spikes_dirs = []
    spikes_bkg_dirs = []
    pattern_dirs = []
    electrodes_dirs = []
    # Format params
    exp, patterns, amplitudes, mice = format_params(exp, patterns, amplitudes, mice) 

    # Radius
    radius = 400

    # Iterate over all parameters and get paths for each combination
    for pattern in patterns:
        for amplitude in amplitudes:
            for mouse in mice:
                root = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

                nodes_dirs.append(os.path.join(root, 'networks', mouse, 'small_network_nodes.h5'))
                spikes_dirs.append(os.path.join(root,  exp, 'output', pattern, amplitude, mouse, 'spikes.csv'))
                #spikes_bkg_dirs.append(os.path.join(root,  exp, 'output', 'bkg', mouse, 'spikes.csv'))
                pattern_dirs.append(os.path.join(root, 'components', 'stimulation', 'patterns', exp, pattern+'.csv'))
                electrodes_dirs.append(os.path.join(root, 'components', 'stimulation', 'electrodes', exp, 'electrodes.csv'))

    return dict(nodes_dirs=nodes_dirs, spikes_dirs=spikes_dirs,spikes_bkg_dirs=spikes_bkg_dirs, pattern_dirs=pattern_dirs, electrodes_dirs=electrodes_dirs, radius=radius, v1=True)

"""
def print_args(names=False, **kwargs):
    print(kwargs)
    length = 1
    for key, value in kwargs.items():
        length *= len(value)
    res = ['']*length
    period = length
    for key, value in kwargs.items():
        period /= len(value)
        iters = length/period
        for i in range(length):
            value_i = (i // period) % len(value)
            value_i = int(value_i)
            if names:
                res[i] += key + str(value[value_i]) + '_'
            else:
                res[i] += str(value[value_i]) + '_'
    for i in range(length):
        res[i] = res[i][:-1]

    print(res)
"""

if __name__ == '__main__':
    # [print(key, ': ', value) for key, value in get_dirs(1,1,10,[0]).items()]
    create_configs(template = r'v1_Anke/exp_1/template.json',
                exp = 'test',
                pattern = 'test', 
                amplitude = 10,
                mice = ['0'], 
                overwrite=True
    )