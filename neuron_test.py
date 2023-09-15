from neuron import h
from neuron.units import ms, mV, µm
h.load_file("stdrun.hoc")

class BallAndStick:
    def __init__(self):
        self._gid = gid
        self.soma = h.Section(name="soma", cell=self)
        self.dend = h.Section(name="dend", cell=self)
        self.dend.connect(self.soma)

    def __repr__(self):
        return "BallAndStick[{}]".format(self._gid)
    
h.topology()
my_cell = BallAndStick(0)
my_other_cell = BallAndStick(1)
h.topology()