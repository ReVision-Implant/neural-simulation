### Imports

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

class CreateWaveform:
    """
    Class to create a voltage waveform and write it to a .csv file.

    This file can be used for BioNet's xstim and comsol modules.
    """

    def __init__(self, piecewise, max_amplitude=None, dt=0.025, path=None, plot=False):
        """ Main function when calling CreateVoltageWaveform(). It 

        :param piecewise: piecewise description of waveform. Each row defines a piece [t_stop, lambda]. t_start of the piece is 0 or t_stop of the previous piece, the lambda expression defines the function. 
        :type piecewise: ndarray
        :param max_amplitude: If specified, normalise waveform to [0,max_amplitude]. Defaults to None.
        :type max_amplitude: int or None, optional
        :param dt: timestep in ms, defaults to 0.025
        :type dt: float, optional
        :param path: if not None, path to the file where the waveform values are saved. Defaults to None.
        :type path: str or None, optional
        :param plot: if true, the waveform shape is plotted. Defaults to False.
        :type plot: bool, optional
        """
    	
        # Initialisation
        self.piecewise = piecewise
        self.max = max_amplitude
        self.dt = dt
        self.path = path
        self.plot = plot
            
        self.create_waveform()

        self.normalise(self.max) if self.max is not None else None

        self.write_to_csv() if self.path is not None else None

        self.plot_waveform() if self.plot is True else None

        return

    def create_waveform(self):
        """ Stores list of timepoints in self.times and calculates corresponding amplitudes, storing them in self.amplitudes. """ 

        self.t_start = 0
        self.times = []
        self.amplitudes = []

        for piece in self.piecewise:
            self.t_stop, self.func = piece
            for t in np.arange(self.t_start, self.t_stop, self.dt):
                self.times.append(t)
                self.amplitudes.append(self.func(t))
            self.t_start = self.t_stop

        self.times = np.round(self.times,10)

    def normalise(self, max=None):
        """ Scale the waveform such that the absolute value of the largest peak is equal to self.max """          

        self.max = max if max is not None else self.max
        assert self.max is not None

        self.amplitudes = self.amplitudes/np.max(np.abs(self.amplitudes))*self.max

    def write_to_csv(self, path=None):
        """ Write the self.times and self.time amplitude to the .csv specified in self.path. """             
        
        self.path = path if path is not None else self.path
        assert os.path.exists(os.path.dirname(self.path))
        
        df = pd.DataFrame({'time': self.times, 'amplitude': self.amplitudes})
        df.to_csv(self.path, sep='\t', index=False)
        return

    def plot_waveform(self):
        """ Plot the waveform stored in self.times and self.amplitudes. """             

        plt.plot(np.append(0,self.times), np.append(0,self.amplitudes)) # Have the line start in (0,0) for visualisation purposes
        plt.xlabel('time [ms]')                         
        plt.ylabel('amplitude [µA]')
        plt.tight_layout()
        plt.show()

        return


def CreateBlockWaveform(n_pulses, phase_1_expr, amp_1_expr, T_1_expr, phase_2_expr, amp_2_expr, T_2_expr):

    piecewise = np.zeros((0,2))
    t_current = 0 
    for i in range(n_pulses):
        
        phase_1 = phase_1_expr(i)
        amp_1 =  amp_1_expr(i)
        T_1 =  T_1_expr(i)
        phase_2 = phase_2_expr(i)
        amp_2 = amp_2_expr(i)
        T_2 = T_2_expr(i)
        
        piecewise_temp1 = [t_current+phase_1, lambda t:amp_1]
        piecewise_temp2 = [t_current+phase_1+T_1, lambda t:0]
        piecewise_temp3 = [t_current+phase_1+T_1+phase_2, lambda t:amp_2]
        piecewise_temp4 = [t_current+phase_1+T_1+phase_2+T_2, lambda t:0]
        
        piecewise = np.vstack((piecewise, piecewise_temp1, piecewise_temp2, piecewise_temp3, piecewise_temp4))

        t_current = t_current+phase_1+T_1+phase_2+T_2

    dir_path = os.path.dirname(os.path.realpath(__file__))
    name = "test.csv"
    path = dir_path + r'/stimulations/' + name
    CreateWaveform(piecewise, max_amplitude = 1, path=path , plot=True)

    return piecewise


if __name__ == '__main__':
    '''
    If you run the file voltage_waveform.py instead of calling if from another file, this part will run.
    
    '''
    piecewise = CreateBlockWaveform(10,
                                phase_1_expr=lambda n:0.1+n/20,
                                amp_1_expr= lambda n:10,
                                T_1_expr=lambda n:0.1,
                                phase_2_expr=lambda n:0.1+n/10,
                                amp_2_expr=lambda n:-5,
                                T_2_expr=lambda n:4.7-3*n/20)