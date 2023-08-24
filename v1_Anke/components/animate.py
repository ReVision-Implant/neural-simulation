import pandas as pd
from hdf5 import HDF5
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.animation import FuncAnimation
import mpl_toolkits.axes_grid1
import matplotlib.widgets
from helper import get_params

class Player(FuncAnimation):
    def __init__(self, fig, func, frames=None, init_func=None, fargs=None,
                 save_count=None, mint=0, maxt=100, pos=(0.125, 0.92), interval=5, **kwargs):
        self.t = 0
        self.min=mint
        self.max=maxt
        self.interval = interval
        self.runs = True
        self.forwards = True
        self.fig = fig
        self.func = func
        self.setup(pos)
        FuncAnimation.__init__(self,self.fig, self.update, frames=self.play(), 
                                           init_func=init_func, fargs=fargs,
                                           save_count=save_count, interval=int(20*interval), **kwargs, )    

    def play(self):
        while self.runs:
            self.t = self.t+self.forwards-(not self.forwards)
            if self.t > self.min and self.t < self.max:
                yield self.t
            else:
                self.stop()
                yield self.t

    def start(self):
        self.runs=True
        self.event_source.start()

    def stop(self, event=None):
        self.runs = False
        self.event_source.stop()

    def forward(self, event=None):
        self.forwards = True
        self.start()
    def backward(self, event=None):
        self.forwards = False
        self.start()
    def oneforward(self, event=None):
        self.forwards = True
        self.onestep()
    def onebackward(self, event=None):
        self.forwards = False
        self.onestep()

    def onestep(self):
        if self.t > self.min and self.t < self.max:
            self.t = self.t+self.forwards-(not self.forwards)
        elif self.t == self.min and self.forwards:
            self.t+=self.interval
        elif self.t == self.max and not self.forwards:
            self.t-=self.interval
        self.func(self.t)
        self.slider.set_val(self.t)
        self.fig.canvas.draw_idle()

    def setup(self, pos):
        playerax = self.fig.add_axes([pos[0],pos[1], 0.64, 0.04])
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(playerax)
        bax = divider.append_axes("right", size="80%", pad=0.05)
        sax = divider.append_axes("right", size="80%", pad=0.05)
        fax = divider.append_axes("right", size="80%", pad=0.05)
        ofax = divider.append_axes("right", size="100%", pad=0.05)
        sliderax = divider.append_axes("right", size="500%", pad=0.07)
        self.button_oneback = matplotlib.widgets.Button(playerax, label='$\u29CF$')
        self.button_back = matplotlib.widgets.Button(bax, label='$\u25C0$')
        self.button_stop = matplotlib.widgets.Button(sax, label='$\u25A0$')
        self.button_forward = matplotlib.widgets.Button(fax, label='$\u25B6$')
        self.button_oneforward = matplotlib.widgets.Button(ofax, label='$\u29D0$')
        self.button_oneback.on_clicked(self.onebackward)
        self.button_back.on_clicked(self.backward)
        self.button_stop.on_clicked(self.stop)
        self.button_forward.on_clicked(self.forward)
        self.button_oneforward.on_clicked(self.oneforward)
        self.slider = matplotlib.widgets.Slider(sliderax, '', 
                                                self.min, 5*self.max, valinit=self.t)
        self.slider.on_changed(self.set_pos)

    def set_pos(self,t):
        self.t = self.slider.val/5
        self.func(self.t)

    def update(self,t):
        self.slider.set_val(5*t)

fig, ax = plt.subplots(figsize=(6,6))
scat = ax.scatter([], [], c="b", s=5,)
ax.set_aspect('equal','box')
ax.set(xlim=[-200,200], ylim=[-200,200], xlabel='x [um]', ylabel='y [um]')
ax.legend()

dictionary = get_params(3, 1, '-', 30, 0)
nodes_dir = dictionary["nodes_dirs"][0]
spikes_dir = dictionary["spikes_dirs"][0]

v1 = True    
if v1 == True:
    node_pos = HDF5(nodes_dir).get_positions_v1()
elif v1 == False:
    node_pos = HDF5(nodes_dir).get_positions()

spikes = pd.read_csv(spikes_dir, sep='\s+')
timestamps = spikes["timestamps"]
node_ids = spikes["node_ids"]
positions = node_pos[node_ids]

circle = positions[:,0]**2 + positions[:,2]**2 <200**2
positions = positions[circle]
timestamps = timestamps[circle]

interval=5
t_max = 100
duration = 10000
frames = int(duration/interval)

def update(t):
    # for each frame, update the data stored on each artist.
    # a = (timestamps<=t)*(timestamps>(t-t_max/frames))
    a = (timestamps<=t*interval)
    x = positions[a,0]
    z = positions[a,2]
    # update the scatter plot:
    data = np.stack([x, z]).T
    scat.set_offsets(data)
    return scat

ani = Player(fig, update, mint=0, maxt=20, interval=interval)

plt.show()