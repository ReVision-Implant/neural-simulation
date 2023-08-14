import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

plt.rcParams.update({'font.size':20})
plt.rcParams.update({'legend.fontsize':30})
plt.rcParams.update({'axes.labelsize':30})
# plt.rcParams.update({'axes.labelweight':'bold'})
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']
hfont = {'fontname':'Verdana'}
plt.rcParams['font.family'] = 'Times New Roman'
mfont = {'fontname':'Times New Roman'}


class Plausibility(): 
        
    def __init__(self, path):

        self.path = path

        # extract useful information from header row in COMSOL output .txt file. 
        header = pd.read_csv(path, sep="\s{3,}", header=None, skiprows=8, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
        header[0] = header[0][2:]                               # remove '% ' before first column name
        if header[3][3] == 'V':
            self._unit = 1000
        elif header[3][3] == 'm':
            self._unit = 1
        for i,col in enumerate(header):                         # remove superfluous characters before actual time value
            if col[0] == "V":
                    header[i] = "V"

        # load data in COMSOL output .txt file.  
        self.comsol = pd.read_csv(path, sep="\s+", header=None, skiprows=9, names=header)           # load data from .txt file

    def plot(self):
        self.data =self.comsol[["z","V"]].sort_values(by=["z"])
        plt.figure(figsize=(8,6))
        z = self.data["z"]-182
        plt.plot(z,-self.data["V"]*10**3)
        plt.plot(z,1/(4*np.pi/3*z)*10**3)
        plt.ylim([0,50])      
        plt.xlim([0,300])
        plt.xticks(np.arange(0,845-182,100))
        plt.yticks(np.arange(0,60,10))
        plt.xlabel('distance from electrode (μm)')
        plt.ylabel('electric potential (mV)')
        plt.legend(['FEM solution','point source'])
        plt.grid(color='lightgrey')
        plt.tight_layout()
        plt.savefig('fig/plausibility.svg')
        plt.savefig('fig/plausibility.png')

Plausibility('plau/plausibility.txt').plot()

class Convergence(): 
        
    def __init__(self, path, type):

        self.path = path
        self.type = type
        
        # extract useful information from header row in COMSOL output .txt file. 
        header = pd.read_csv(path, sep="\s{1,}", header=None, skiprows=5, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
        header[0] = 0                       # remove '% ' before first column name
        self.header = header
        self.n = int(len(self.header)/2)
        for i in range(1,len(self.header)):
            self.header[i] = int(i)
        self.legend = [1,2,3,4,5,6,7,8]
        if self.type == 'g':
            self.flag = 0*self.n
        elif self.type == '-':
            self.flag = 1*self.n-1
        # load data in COMSOL output .txt file.  
        self.comsol = pd.read_csv(path, sep="\s+", header=None, skiprows=5, names=self.header)           # load data from .txt file
        self.dofs = np.array([3.4982*10**5, 1.0174*10**6, 1.9807*10**6, 2.8482*10**6, 3.3737*10**6])/10**6
        

    def plotV_max(self, ax):

        for i in range(1+self.flag,self.n+1+self.flag):
            self.data =self.comsol[[0,int(i)]].sort_values(by=[0])
            ax.plot(self.dofs,self.data[int(i)]*10**3)
        ax.set_xscale('linear')
        # ax.set_xlabel('number of degrees of freedom ($10^6$)')
        ax.set_ylabel('maximum potential (mV)')
        # ax.title('Convergence of the maximum potential')
        # ax.savefig('V_max' + self.type + '.svg')

    def plotA(self, ax):

        self.comsol[[7,8]] = self.comsol[[33,34]] 
        self.comsol = self.comsol[np.arange(33)]
        for i in range(1+self.flag,self.n+self.flag,2):
            self.data =self.comsol[[0,int(i)]].sort_values(by=[0])
            ax[0].plot(self.dofs,self.data[int(i)]*10**12)
        ax[0].set_xscale('linear')
        # ax[0].set_xlabel('number of degrees of freedom ($10^6$)')
        ax[0].set_ylabel('∫ V>1mV (10⁶ μm³)')
        # ax.title('Convergence of the maximum potential')
        # ax[0].savefig('AV' + self.type + '.svg')
        
        for i in range(2+self.flag,self.n+self.flag,2):
            self.data = self.comsol[[0,int(i)]].sort_values(by=[0])
            ax[1].plot(self.dofs,self.data[int(i)]*10**12)
        ax[1].set_xscale('linear')
        # ax[1].set_xlabel('number of degrees of freedom ($10^6$)')
        ax[1].set_ylabel('∫ J>0.05A/m² (10⁶ μm³)')
        # ax.title('Convergence of the maximum potential')
        # ax[1].savefig('AJ' + self.type + '.svg')

    def plotI(self, ax):

        for i in range(1,2*self.n+1):
            self.data =self.comsol[[0,int(i)]].sort_values(by=[0])
            ax.plot(self.dofs,self.data[int(i)]*10**9)
        ax.set_xscale('linear')
        ax.set_xlabel('number of degrees of freedom (10⁶)')
        ax.set_ylabel('return current (pA)')
        # ax.title('Convergence of the maximum potential')
        # ax.savefig('I' + self.type + '.svg')

    
    def plotV(self, ax):

        for i in range(1,2*self.n+1):
            self.data =self.comsol[[0,int(i)]].sort_values(by=[0])
            ax.plot(self.dofs,self.data[int(i)]*10**3)
        ax.set_xscale('linear')
        ax.set_xlabel('number of degrees of freedom (10⁶)')
        ax.set_ylabel('return potential (mV)')
        # ax.title('Convergence of the maximum potential')
        # ax.savefig('V' + self.type + '.svg')



fig, axs = plt.subplots(4,2, sharex=True, figsize = (16,20))
axs[0,0].set_ylim([46.25,48.75])
axs[0,1].set_ylim([46.25,48.75])
axs[3,1].set_ylim([-48.75,-46.25])
Convergence('mesh/V_max.txt','g').plotV_max(axs[0,0])
Convergence('mesh/V_max.txt','-').plotV_max(axs[0,1])
Convergence('mesh/A.txt','g').plotA(axs[[1,2],0])
Convergence('mesh/A.txt','-').plotA(axs[[1,2],1])
Convergence('mesh/I.txt','g').plotI(axs[3,0])
Convergence('mesh/V.txt','-').plotV(axs[3,1])
for ax in axs.flatten():
    ax.get_yaxis().set_label_coords(-0.15,0.5)
    ax.get_xaxis().set_label_coords(0.5,-0.2)
    ax.axvline(x=1.0174, c='k', ls='--', lw=0.5)
    ax.grid(True, axis='y', color='lightgrey')
axs[3,0].set_xticks(np.round(np.array([3.4982*10**5, 1.0174*10**6, 1.9807*10**6, 2.8482*10**6, 3.3737*10**6])/10**6,2))
axs[0,0].set_yticks(np.arange(46.5,49,0.5))
axs[0,1].set_yticks(np.arange(46.5,49,0.5))
axs[3,1].set_yticks(np.arange(-48.5,-46,0.5))
axs[0,0].set_title('return electrode: ground', y=1.1, size=35, weight='bold')
axs[0,1].set_title('return electrode: current', y=1.1, size=35, weight='bold')
fig.suptitle(' ', size=40)
fig.legend(range(1,9), loc='upper center', ncol=8, bbox_to_anchor=(0.5, 1))
fig.tight_layout(w_pad=2, h_pad=2)
plt.savefig('fig/mesh.svg')
plt.savefig('fig/mesh.png')
# plt.show()


class External(): 
        
    def __init__(self, path):

        self.path = path
        self.type = type
        
        # extract useful information from header row in COMSOL output .txt file. 
        header = pd.read_csv(path, sep="\s{1,}", header=None, skiprows=5, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
        header[0] = 0                       # remove '% ' before first column name
        self.header = header
        self.n = int(len(self.header)/2)
        for i in range(1,len(self.header)):
            self.header[i] = int(i)
        self.legend = ['ground', 'current']

        # load data in COMSOL output .txt file.  
        self.comsol = pd.read_csv(path, sep="\s+", header=None, skiprows=5, names=self.header)           # load data from .txt file
        self.dofs = np.array([1,2,5,10,20,50])
        self.comsol = self.comsol.iloc[:-2,:]

    def plotV_max(self, ax):

        ax.plot(self.dofs,self.comsol[1]*10**3)
        ax.plot(self.dofs,self.comsol[2]*10**3)
        ax.set_xscale('log')
        ax.set_ylabel('maximum potential (mV)')

    def plotA(self, ax):

        ax[0].plot(self.dofs,self.comsol[1]*10**12)
        ax[0].plot(self.dofs,self.comsol[3]*10**12)
        ax[0].set_xscale('log')
        ax[0].set_ylabel('∫ V>1mV (10⁶ μm³)')

        ax[1].plot(self.dofs,self.comsol[2]*10**12)
        ax[1].plot(self.dofs,self.comsol[4]*10**12)
        ax[1].set_xscale('log')
        ax[1].set_ylabel('∫ J>0.05A/m² (10⁶ μm³)')

    def plotIV(self, ax):

        ax.plot(self.dofs,self.comsol[1]*10**9)
        ax.set_xscale('log')
        ax.set_ylabel('return current (pA)', color='#1f77b4')
        ax2 = ax.twinx()
        ax2.plot(self.dofs,self.comsol[2]*10**3, color='tab:orange')
        ax2.set_ylabel('return voltage (mV)', color='tab:orange', labelpad=20)
        ax2.set_yticks(np.arange(-47.3,-47.0,0.1))        


fig, axs = plt.subplots(2,2, sharex=True, figsize = (16,10))
External('ext/V_max.txt').plotV_max(axs[0,0])
External('ext/A.txt').plotA([axs[0,1],axs[1,0]])
External('ext/IV.txt').plotIV(axs[1,1])
for ax in axs.flatten():
    ax.get_yaxis().set_label_coords(-0.15,0.5)
    ax.get_xaxis().set_label_coords(0.5,-0.2)
    ax.axvline(x=2, c='k', ls='--', lw=0.5)
    ax.grid(True, axis='y', color='lightgrey')
fig.suptitle(' ', size=40)
fig.legend(['ground', 'current'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1))
fig.tight_layout(w_pad=2, h_pad=2)
plt.savefig('fig/ext.png')
plt.savefig('fig/ext.svg')
# plt.show()


class Stability(): 
        
    def __init__(self, path):

        self.path = path
        self.type = type
        
        # extract useful information from header row in COMSOL output .txt file. 
        header = pd.read_csv(path, sep="\s{1,}", header=None, skiprows=5, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
        header[0] = 0                       # remove '% ' before first column name
        self.header = header
        self.n = int(len(self.header)/2)
        for i in range(1,len(self.header)):
            self.header[i] = int(i)
        self.legend = ['ground', 'current']

        # load data in COMSOL output .txt file.  
        self.comsol = pd.read_csv(path, sep="\s+", header=None, skiprows=5, names=self.header)           # load data from .txt file
        self.IrOx = np.array([10**i for i in range(4,11)])
        self.PI = np.array([10**i for i in range(-10,-3,1)])

    def plotV_max(self, ax):

        ax[0].plot(self.IrOx,self.comsol[1]*10**3)
        ax[0].plot(self.IrOx,self.comsol[2]*10**3)
        ax[0].set_xscale('log')
        ax[0].set_ylabel('maximum potential (mV)')

        ax[1].plot(self.PI,self.comsol[3]*10**3)
        ax[1].plot(self.PI,self.comsol[4]*10**3)
        ax[1].set_xscale('log')
        ax[1].set_ylabel('maximum potential (mV)')

    def plotA(self, ax):

        ax[0,0].plot(self.IrOx,self.comsol[1]*10**12)
        ax[0,0].plot(self.IrOx,self.comsol[3]*10**12)
        ax[0,0].set_xscale('log')
        ax[0,0].set_ylabel('∫ V>1mV (10⁶ μm³)')

        ax[0,1].plot(self.PI,self.comsol[5]*10**12)
        ax[0,1].plot(self.PI,self.comsol[7]*10**12)
        ax[0,1].set_xscale('log')
        ax[0,1].set_ylabel('∫ V>1mV (10⁶ μm³)')

        ax[1,0].plot(self.IrOx,self.comsol[2]*10**12)
        ax[1,0].plot(self.IrOx,self.comsol[4]*10**12)
        ax[1,0].set_xscale('log')
        ax[1,0].set_ylabel('∫ J>0.05A/m² (10⁶ μm³)')

        ax[1,1].plot(self.PI,self.comsol[6]*10**12)
        ax[1,1].plot(self.PI,self.comsol[8]*10**12)
        ax[1,1].set_xscale('log')
        ax[1,1].set_ylabel('∫ J>0.05A/m² (10⁶ μm³)')
        
    def plotIV(self, ax):

        ax[0].plot(self.IrOx,self.comsol[1]*10**9)
        ax[0].set_xscale('log')
        ax[0].set_ylabel('return current (pA)', color='#1f77b4')
        ax[0].set_yticks(np.arange(-16,-15.75,0.1))        
        ax2 = ax[0].twinx()
        ax2.plot(self.IrOx,self.comsol[2]*10**3, color='tab:orange')
        ax2.set_ylabel('return voltage (mV)', color='tab:orange', labelpad=20)
        ax2.set_yticks(np.arange(-47.4,-47.1,0.1))        

        ax[1].plot(self.PI,self.comsol[3]*10**9)
        ax[1].set_xscale('log')
        ax[1].set_ylabel('return current (pA)', color='#1f77b4')   
        ax[1].set_yticks(np.arange(-16,-15.75,0.1))           
        ax1 = ax[1].twinx()
        ax1.plot(self.PI,self.comsol[2]*10**3, color='tab:orange')
        ax1.set_ylabel('return voltage (mV)', color='tab:orange', labelpad=20)
        ax1.set_yticks(np.arange(-47.4,-47.1,0.1))        

        ax[0].set_xlabel(r'$\sigma_{IrOx}$')
        ax[1].set_xlabel(r'$\sigma_{PI}$')

fig, axs = plt.subplots(4,2, figsize = (16,20))
Stability('stab/V_max.txt').plotV_max(axs[0,:])
Stability('stab/A.txt').plotA(axs[[1,2],:])
Stability('stab/IV.txt').plotIV(axs[3,:])
for ax in axs.flatten():
    ax.get_yaxis().set_label_coords(-0.2,0.5)
    ax.get_xaxis().set_label_coords(0.5,-0.1)
    ax.grid(True, axis='y', color='lightgrey')
for ax in axs[:3,:].flatten():
    ax.set_xticklabels([])
for ax in axs[:,0].flatten():
    ax.axvline(x=2.5e6, c='k', ls='--', lw=0.5)
for ax in axs[:,1].flatten():
    ax.axvline(x=1e-7, c='k', ls='--', lw=0.5)

axs[0,0].set_title('IrOx', y=1.2, size=35, weight='bold')
axs[0,1].set_title('PI', y=1.2, size=35, weight='bold')
fig.suptitle(' ', size=40)
fig.legend(['ground', 'current'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1))
fig.tight_layout(w_pad=2, h_pad=2)
plt.savefig('fig/stab.png')
plt.savefig('fig/stab.svg')
# plt.show()


class Validation(): 
        
    def __init__(self, path):

        self.path = path

        # extract useful information from header row in COMSOL output .txt file. 
        header = pd.read_csv(path, sep="\s{3,}", header=None, skiprows=8, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
        header[0] = header[0][2:]                               # remove '% ' before first column name
        if header[3][3] == 'V':
            self._unit = 1000
        elif header[3][3] == 'm':
            self._unit = 1
        for i,col in enumerate(header):                         # remove superfluous characters before actual time value
            if col[0] == "V":
                    header[i] = "V"

        # load data in COMSOL output .txt file.  
        self.comsol = pd.read_csv(path, sep="\s+", header=None, skiprows=9, names=header)           # load data from .txt file

    def plot(self):
        self.data =self.comsol[["z","V"]].sort_values(by=["z"])
        plt.figure(figsize=(8,6))
        z = self.data["z"]-182
        plt.plot(z,-self.data["V"]*10**3/10)
        plt.plot(z,1/(4*np.pi/3*z)*10**3/10)
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim([0.01,100])      
        plt.xlim([1,800])
        # plt.xticks(np.arange(0,200,25))
        # plt.yticks(np.arange(0,60,10))
        plt.xlabel('distance from electrode (μm)')
        plt.ylabel('electric potential (mV)')
        plt.legend(['FEM solution','point source'])
        plt.grid(color='lightgrey')
        plt.tight_layout()
        plt.savefig('fig/val.svg')
        plt.savefig('fig/val.png')

Validation('plau/plausibility.txt').plot()
# plt.show()