import json
import os

def create_configs(template, exp, electrodes, stim_types, amplitudes, networks, overwrite=False):
    
    
    # Load template config.json
    with open(template, 'r') as read_file:
        data = json.load(read_file)

    # Nested loop over all specified parameters (converted to str)
    for electrode in [str(i) for i in electrodes]:
        for stim_type in [str(i) for i in stim_types]:

            # Get path to comsol.txt and set it in config.json
            comsol_file = '$STIM_DIR/' + exp + '/0' + electrode + stim_type + '.txt'
            data['inputs']['Extracellular_Stim']['comsol_file'] = comsol_file

            for amplitude in [str(i) for i in amplitudes]:
                
                # Set amplitude in config.json
                data['inputs']['Extracellular_Stim']['amplitude'] = int(amplitude)

                for network in [str(i) for i in networks]:
                    
                    # Set network dir in config.json
                    data['manifest']['$NETWORK_DIR'] = data['manifest']['$NETWORK_DIR'][0:-1] + network
                                                
                    # Get output dir name (same as config.json name) and set it in config.json
                    name = '0' + electrode + stim_type + '_' + amplitude + '_' + network
                    data['manifest']["$OUTPUT_DIR"] = '$BASE_DIR/' + exp + '/output/' + electrode + '/' + stim_type + '/' + amplitude + '/' + name

                    # Get path for config.json and create directory if it does not exist
                    root = os.path.dirname(os.path.abspath(__file__))
                    print(root)
                    path = root + '/' + exp + '/config/' + electrode + '/' + stim_type + '/' + amplitude
                    if not os.path.exists(path):
                        os.makedirs(path)
                    
                    # Save config.json unless it already exists and overwrite is False     
                    if (overwrite == True) or (not os.path.exists(path + '/config_' + name + '.json')):
                        with open(path + '/config_' + name + '.json', 'w') as write_file:
                            json.dump(data, write_file, indent=2)
                            print(name)
    return data

def delete_configs(exp, electrodes, stim_types, amplitudes, networks):

    # Nested loop over all specified parameters (converted to strings)
    for electrode in [str(i) for i in electrodes]:

        for stim_type in [str(i) for i in stim_types]:

            for amplitude in [str(i) for i in amplitudes]:

                for network in [str(i) for i in networks]:
                    
                    # Get file names and paths corresponding to specified parameters
                    root = os.path.dirname(os.path.abspath(__file__))
                    name = '0' + electrode + stim_type + '_' + amplitude + '_' + network
                    path = root + '/' + '/config/'+ electrode + '/' + stim_type + '/' + amplitude

                    # Remove file if it exists
                    if os.path.exists(path + '/config_' + name + '.json'):
                        os.remove(path + '/config_' + name + '.json')
                        print('removed file: ', name)

    # Remove empty directories
    walk = list(os.walk('./config'))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)
            print('removed folder: ', path)

