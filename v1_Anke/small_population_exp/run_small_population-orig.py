#choose correct folder depending on directed axons or just long axons
#make sure config file also includes a csv file

from bmtk.simulator import bionet

conf = bionet.Config.from_json('simulation/config.json')
#conf=bionet.Config.from_json('simulation_long_axons/config.json')
#conf=bionet.Config.from_json(simulation_directed_axons/config.json')
conf.build_env()
net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)
sim.run()