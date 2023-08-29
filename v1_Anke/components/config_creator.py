import json
import os
import pandas as pd

def create_configs(template, exp, patterns, amplitudes, networks, overwrite=False):

    # Convert integer/string to list if necessary
    patterns = [patterns] if not isinstance(patterns, list) else patterns
    amplitudes = [amplitudes] if not isinstance(amplitudes, list) else amplitudes
    networks = [networks] if not isinstance(networks, list) else networks
    exp = 'exp' + str(exp) if str(exp)[:3] != 'exp' else str(exp)

    # Load template config.json
    with open(template, 'r') as read_file:
        data = json.load(read_file)

    # Nested loop over all specified parameters (converted to str)
    for pattern in ['pattern'+str(i) for i in patterns]:

        # Get pattern parameters from pattern.csv
        stim_params = get_stim_params(exp, pattern)

        # Set paths to comsol.txt and waveform.csv files in config.json
        comsol_files = [os.path.join('$STIM_DIR', 'comsol', file) for file in stim_params['comsol_files']]
        waveforms = [os.path.join('$STIM_DIR', 'waveforms', file) for file in stim_params['waveforms']]
        data['inputs']['Extracellular_Stim']['comsol_files'] = comsol_files
        data['inputs']['Extracellular_Stim']['waveforms'] = waveforms

        for amplitude in [i for i in amplitudes]:
            
            # Set amplitude in config.json
            amplitudes = [i*float(amplitude) for i in stim_params['amplitudes']]
            data['inputs']['Extracellular_Stim']['amplitudes'] = amplitudes
            amplitude = str(amplitude)+'uA'

            for network in ['network'+str(i) for i in networks]:
                
                # Set network dir in config.json
                network_dir = os.path.join(os.path.dirname(data['manifest']['$NETWORK_DIR']), network)
                data['manifest']['$NETWORK_DIR'] = network_dir
                                            
                # Set output dir in config.json
                output_dir = os.path.join('$BASE_DIR', exp , 'output', pattern, amplitude, network)
                data['manifest']["$OUTPUT_DIR"] = output_dir

                # Get path for config.json and create directory if it does not exist
                root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                file_dir = os.path.join(root, exp, 'config', pattern, amplitude)
                if not os.path.exists(file_dir):
                    os.makedirs(file_dir)
                
                # Save config.json unless it already exists and overwrite is False
                file_name =  'config.' + exp[0]+exp[-1] + '_' + pattern[0]+pattern[-1] + '_' + amplitude[0]+amplitude[-1] + '_' + network[0]+network[-1] + '.json' 
                print(os.path.join(file_dir, file_name))
                if (overwrite == True) or (not os.path.exists(os.path.join(file_dir, file_name))):
                    with open(os.path.join(file_dir, file_name), 'w') as write_file:
                        json.dump(data, write_file, indent=2)
                        print('created file: ',file_name)
    return data

def delete_configs(exp, patterns=None, amplitudes=None, networks=None):

    # Convert integer/string to list if necessary
    patterns = [patterns] if not isinstance(patterns, list) else patterns
    amplitudes = [amplitudes] if not isinstance(amplitudes, list) else amplitudes
    networks = [networks] if not isinstance(networks, list) else networks
    exp = 'exp' + str(exp) if str(exp)[:3] != 'exp' else str(exp)

    # 
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 

    # If any of the input arguments is None, skip this statement and just remove empty directories
    if patterns is None or amplitudes is None or networks is None:
        pass
    else:
        # Nested loop over all specified parameters (converted to strings)
        for pattern in ['pattern'+str(i) for i in patterns]:

            for amplitude in [str(i)+'uA' for i in amplitudes]:

                for network in ['network'+str(i) for i in networks]:
                    
                    # Get file name and path 
                    file_path = os.path.join(root, exp, 'config', pattern, amplitude)                
                    file_name = 'config.' + exp[0]+exp[-1] + '_' + pattern[0]+pattern[-1] + '_' + amplitude[0]+amplitude[-1] + '_' + network[0]+network[-1] + '.json'
                    
                    # Remove file if it exists
                    if os.path.exists(os.path.join(file_path, file_name)):
                        os.remove(os.path.join(file_path, file_name))
                        print('removed file: ', file_name)

    # Remove empty directories
    walk = list(os.walk(os.path.join(root,exp,'config')))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)
            print('removed folder: ', path)

def get_stim_params(exp, pattern):

    exp = 'exp'+str(exp) if str(exp)[:3] != 'exp' else str(exp)
    pattern = str(pattern)+'.csv' if str(pattern)[-4:] != '.csv' else pattern
    pattern = 'pattern'+pattern if pattern[:7] != 'pattern' else pattern
    root = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(root, 'stimulation', 'patterns', exp, pattern)
    df = pd.read_csv(path, sep='\s+')
    
    comsol_files = df['comsol'].tolist()
    waveforms = df['waveform'].tolist()
    amplitudes = df['amplitude'].tolist()
    
    return dict(comsol_files=comsol_files, waveforms=waveforms, amplitudes=amplitudes)

amplitudes = get_stim_params(1, 1)['amplitudes']
# print(amplitudes*10)

create_configs(r'..\v1_Anke\exp1\template.json', 'exp1', 1, 10, 0, overwrite=True)
# delete_configs('exp1', 1, 10, 0)