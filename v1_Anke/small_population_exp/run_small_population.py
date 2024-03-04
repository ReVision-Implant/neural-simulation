from bmtk.simulator import bionet


conf = bionet.Config.from_json('simulation/config.json')
conf.build_env()
net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)
sim.run()