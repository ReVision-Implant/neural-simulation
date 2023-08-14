import json
import os

with open('config/1/-/10/config_01-_10_0.json', 'r') as read_file:
    data = json.load(read_file)

def create_configs(electrodes, stim_types, amplitudes, networks, overwrite=False):
    for electrode in electrodes:

        for stim_type in stim_types:

            comsol_file = '$STIM_DIR/exp3/0' + str(electrode) + str(stim_type) + '.txt'
            data['inputs']['Extracellular_Stim']['comsol_file'] = comsol_file

            for amplitude in amplitudes:

                data['inputs']['Extracellular_Stim']['amplitude'] = amplitude

                for network in networks:

                    name = '0' + str(electrode) + str(stim_type) + '_' + str(amplitude) + '_' + str(network)

                    data['manifest']['$NETWORK_DIR'] = data['manifest']['$NETWORK_DIR'][0:-1] + str(network)
                    data['manifest']["$OUTPUT_DIR"] = '$BASE_DIR/exp3/output/' + str(electrode) + '/' + str(stim_type) + '/' + str(amplitude) + '/' + name

                    path = 'config/'+ str(electrode) + '/' + str(stim_type) + '/' + str(amplitude)
                    if not os.path.exists(path):
                        os.makedirs(path)
                    if (overwrite == True) or (not os.path.exists(path + '/config_' + name + '.json')):
                        with open(path + '/config_' + name + '.json', 'w') as write_file:
                            json.dump(data, write_file, indent=2)
                            print(name)


def delete_configs(electrodes, stim_types, amplitudes, networks):
    for electrode in electrodes:

        for stim_type in stim_types:

            for amplitude in amplitudes:

                for network in networks:

                    name = '0' + str(electrode) + str(stim_type) + '_' + str(amplitude) + '_' + str(network)
                    path = 'config/'+ str(electrode) + '/' + str(stim_type) + '/' + str(amplitude)
                    if os.path.exists(path + '/config_' + name + '.json'):
                        os.remove(path + '/config_' + name + '.json')
                        print(name, 'removed')


    walk = list(os.walk('./config'))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)


delete_configs(
    electrodes = range(1,4),
    stim_types = ['-','g'],
    amplitudes = [10,20,30],
    networks = range(0,3),
)

create_configs(
    electrodes = range(1,4),
    stim_types = ['-','g'],
    amplitudes = [10,20,30],
    networks = range(0,3),
)

# create_configs(
#     electrodes = range(1,9),
#     stim_types = ['-'],
#     amplitudes = [10,20,30],
#     networks = range(0,3)
# )


