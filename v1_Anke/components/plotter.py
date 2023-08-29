import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, pearsonr
from plotter_helper import get_image, plot_image
from spikes_helper import get_params, get_spikes

def fractions(FTA, FNA, AR, save=None):

    fig, ax = plt.subplots(figsize=(8,6))
    
    n_columns = np.size(FTA,1)

    FTA = 100*np.vstack((np.zeros((1,n_columns)),FTA))
    error_low = np.mean(FTA,1) - np.min(FTA, 1)
    error_high = np.max(FTA,1) - np.mean(FTA,1)
    amplitudes = [0,10,20,30]
    p = ax.errorbar(amplitudes, np.mean(FTA,1), yerr=np.vstack((error_low, error_high)), capsize=5, label='FTA (%)')

    FNA = 100*np.vstack((np.zeros((1,n_columns)),FNA))
    error_low = np.mean(FNA,1) - np.min(FNA, 1)
    error_high = np.max(FNA,1) - np.mean(FNA,1)
    p2 = ax.errorbar(amplitudes, np.mean(FNA,1), yerr=np.vstack((error_low, error_high)), capsize=5, linestyle='--', c='red', label='FNA (%)')

    ax.set_ylabel('fraction (%)', labelpad=2)
    ax.set_xlabel('amplitude (μA)')
    ax.set_xticks(amplitudes)
    if save[11] == '1':
        ax.set_ylim([0,2])
    elif save[11] == '3':
        ax.set_ylim([0,100])

    # ax2 = ax.twinx()
    AR = np.vstack((np.zeros((1,n_columns)),AR))
    # error_low = np.mean(AR,1) - np.min(AR, 1)
    # error_high = np.max(AR,1) - np.mean(AR,1)
    # p3 = ax2.errorbar(amplitudes, np.mean(AR,1), yerr=np.vstack((error_low, error_high)), capsize=5, linestyle='-.', c='green', label='AR (um)')

    # ax2.set_ylabel('FNA')
    # if save[11] == '1':
    #     ax2.set_ylim([0,200])
    #     ax2.set_ylabel('average radius (μm)')
    # elif save[11] == '3':
    #     ax2.set_ylim([0,200])  
    #     ax2.set_ylabel('average radius (μm)')
    ax.legend(handles=[p,p2], labels=['FTA','FNA'], loc='upper left')
    fig.tight_layout()
    if save is not None:
        fig.savefig(save+'.svg')
        fig.savefig(save+'.png')
    else:
        plt.show()

    r_FTA = pearsonr(np.mean(FTA,1), amplitudes)
    r_FNA = pearsonr(np.mean(FNA,1), amplitudes)
    r_AR = pearsonr(np.mean(AR,1), amplitudes)

    return r_FTA, r_FTA.confidence_interval(confidence_level=0.95), r_FNA, r_FNA.confidence_interval(confidence_level=0.95)

def plot_3x3(exp, stim_type, save=None):
    radius = 400 if exp == 1 else 200
    ### Initialisation
    images = np.zeros((3,3,radius*2,radius*2))
    electrodes = [1,2,3]
    amplitudes = [10,20,30]
    FTA = np.zeros((3,3))
    FNA = np.zeros((3,3))
    AR = np.zeros((3,3))
    fig, axs = plt.subplots(3,3, sharex=True, sharey=True, figsize=(16,13))
    for column,electrode in enumerate(electrodes):
        for row,amplitude in enumerate(amplitudes):
            images[row,column,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
            node_pos = np.zeros((0,3))
            n_spikes = np.zeros(0)
            for network in ['0','1','2']:
                node_pos_temp, n_spikes_temp = get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))
                node_pos = np.vstack((node_pos, node_pos_temp))
                n_spikes = np.hstack((n_spikes, n_spikes_temp))
            FNA[row,column] = np.sum(n_spikes>0)/np.sum(n_spikes>-1)
            radius = np.sqrt(node_pos[:,0]**2 + node_pos[:,2]**2)
            AR[row, column] = np.average(radius, weights= n_spikes)
    # images = np.array(images)/np.nanmax(images)*0.9
    # images = np.log10(images+0.1)+1

    images = np.array(images)/np.nanmax(images)

    for i in range(np.size(images,0)):
        for j in range(np.size(images,1)):
            image = images[i,j,:,:]/np.nanmax(images)
            FTA[i,j] = np.sum(image>0.001)/np.sum(image!=None)    

    for column,electrode in enumerate(electrodes):
        for row,amplitude in enumerate(amplitudes):
            image = images[row,column,:,:]
            im = plot_image(axs[row,column], image, **get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
            if row == 0:
                axs[row,column].set_title('return electrode ' + str(column+1), y=1.1, size=30)
            if column == 0:
                axs[row,column].set_ylabel(str(amplitude) + ' μA', rotation=90, size=30)
                axs[row,column].yaxis.set_label_coords(-.3,.5)
    handles, labels = axs[0,0].get_legend_handles_labels()
    fig.tight_layout(h_pad=0.0,w_pad=.7)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=30, shrink=1)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, -0.1), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()

    return FTA, FNA, AR

def plot_2x3(exps, stim_types, save=None):
    ### Initialisation
    radius = 200 if len(exps)>1 else 400
    images = np.zeros((3,3,radius*2,radius*2))
    amplitude = 20
    electrodes = [1,2,3]
    fig, axs = plt.subplots(2,3, sharex=True, sharey=True, figsize=(16,10))
    for column,electrode in enumerate(electrodes):
        for row1, exp in enumerate(exps):
            for row2, stim_type in enumerate(stim_types):                
                row = row1 + row2
                images[row,column,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
                node_pos = np.zeros((0,3))
                n_spikes = np.zeros(0)
                for network in ['0','1','2']:
                    node_pos_temp, n_spikes_temp = get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))
                    node_pos = np.vstack((node_pos, node_pos_temp))
                    n_spikes = np.hstack((n_spikes, n_spikes_temp))
                radius = np.sqrt(node_pos[:,0]**2 + node_pos[:,2]**2)

    images = np.array(images)/np.nanmax(images)

    for column,electrode in enumerate(electrodes):
        for row1, exp in enumerate(exps):
            for row2, stim_type in enumerate(stim_types):
                row = row1 + row2
                image = images[row,column,:,:]
                im = plot_image(axs[row,column], image, **get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
                if row == 0:
                    axs[row,column].set_title('return electrode ' + str(column+1), y=1.1, size=30)
                if column == 0:
                    if len(exps) > 1:
                        if len(exp) == 1:
                            axs[row,column].set_ylabel('cathodic-first', y=1.1, size=30)              
                        elif len(exp) == 2:
                            axs[row,column].set_ylabel('anodic-first', y=1.1, size=30) 
                    elif len(stim_types) > 1:
                        if stim_type == '-':
                            axs[row,column].set_ylabel('active return', y=1.1, size=30)              
                        elif stim_type == 'g':
                            axs[row,column].set_ylabel('grounded return', y=1.1, size=30)                     
                    axs[row,column].yaxis.set_label_coords(-.3,.5)
    handles, labels = axs[0,0].get_legend_handles_labels()
    fig.tight_layout(h_pad=-3,w_pad=.7)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=27, shrink=0.9)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, -0.1), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()


def plot_2x2(exps, stim_types, amplitude, save=None):
    ### Initialisation
    radius = 200 if len(exps)>1 else 400
    images = np.zeros((3,3,radius*2,radius*2))
    electrodes = [1,2]
    fig, axs = plt.subplots(2,2, sharex=True, sharey=True, figsize=(11,9))
    for column,electrode in enumerate(electrodes):
        for row1, exp in enumerate(exps):
            for row2, stim_type in enumerate(stim_types):                
                row = row1 + row2
                images[row,column,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
                node_pos = np.zeros((0,3))
                n_spikes = np.zeros(0)
                for network in ['0','1','2']:
                    node_pos_temp, n_spikes_temp = get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))
                    node_pos = np.vstack((node_pos, node_pos_temp))
                    n_spikes = np.hstack((n_spikes, n_spikes_temp))
                radius = np.sqrt(node_pos[:,0]**2 + node_pos[:,2]**2)

    images = np.array(images)/np.nanmax(images)

    for column,electrode in enumerate(electrodes):
        for row1, exp in enumerate(exps):
            for row2, stim_type in enumerate(stim_types):
                row = row1 + row2
                image = images[row,column,:,:]
                im = plot_image(axs[row,column], image, **get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
                if row == 0:
                    axs[row,column].set_title('return electrode ' + str(column+1), y=1.1, size=30)
                if column == 0:
                    if len(exps) > 1:
                        if exp[0] == '3':
                            axs[row,column].set_ylabel('original \norientaiton', y=1.1, size=30)              
                        elif exp[0] == '2':
                            axs[row,column].set_ylabel('opposite \norientation', y=1.1, size=30) 
                    elif len(stim_types) > 1:
                        if stim_type == '-':
                            axs[row,column].set_ylabel('active return', y=1.1, size=30)              
                        elif stim_type == 'g':
                            axs[row,column].set_ylabel('grounded return', y=1.1, size=30)                     
                    axs[row,column].yaxis.set_label_coords(-.3,.5)
    handles, labels = axs[0,0].get_legend_handles_labels()
    fig.tight_layout(h_pad=0,w_pad=1.7)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=30, shrink=1)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, -0.1), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()


def plot_app(exp, save=None):
    ### Initialisation
    radius = 200
    images = np.zeros((2,3,radius*2,radius*2))
    electrodes = [1,2]
    amplitudes = [10,20,30]
    stim_type = '-'
    fig, axs = plt.subplots(2,3, sharex=True, sharey=True, figsize=(16,10))
    for row,electrode in enumerate(electrodes):
        for column, amplitude in enumerate(amplitudes):
            images[row,column,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
            node_pos = np.zeros((0,3))
            n_spikes = np.zeros(0)
            for network in ['0','1','2']:
                node_pos_temp, n_spikes_temp = get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))
                node_pos = np.vstack((node_pos, node_pos_temp))
                n_spikes = np.hstack((n_spikes, n_spikes_temp))
            radius = np.sqrt(node_pos[:,0]**2 + node_pos[:,2]**2)

    images = np.array(images)/np.nanmax(images)

    for row,electrode in enumerate(electrodes):
        for column, amplitude in enumerate(amplitudes):
            image = images[row,column,:,:]
            im = plot_image(axs[row,column], image, **get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
            if row == 0:
                axs[row,column].set_title(str(amplitude)+'μA', y=1.1, size=30)
            if column == 0:
                axs[row,column].set_ylabel('return electrode ' + str(row+1), size=30)              
            axs[row,column].yaxis.set_label_coords(-.3,.5)


    handles, labels = axs[0,0].get_legend_handles_labels()
    fig.tight_layout(h_pad=-1.4,w_pad=.4)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=27.9, shrink=0.93)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, -0.1), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()


def plot_1x3(exps, electrodes, stim_types, amplitudes, networks, save=None, labels=None, radius=200):

    ### Initialisation
    images = np.zeros((3,radius*2,radius*2))
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True, figsize=(16,7))

    for i,exp in enumerate(exps):
        for j,electrode in enumerate(electrodes):
            for k,stim_type in enumerate(stim_types):
                for l,amplitude in enumerate(amplitudes):
                    if stim_type in [0,'0']:
                        exp = 0
                        stim_type = 'g'
                        electrode_temp = 0
                    else:
                        electrode_temp = electrode
                        exp = 3
                    if len(networks) == 3:
                        n = i+j+k+l
                        images[n,:,:] = get_image(**get_params(exp, electrode_temp, stim_type, amplitude, ['0','1','2']))
                    elif len(networks) == 2:
                        for m,network in enumerate(networks):
                            images[m,:,:] = get_image(**get_params(exp, electrode_temp, stim_type, amplitude, [network]))

    images = np.array(images)/np.nanmax(images)

    for i,exp in enumerate(exps):
        for j,electrode in enumerate(electrodes):
            for k,stim_type in enumerate(stim_types):
                for l,amplitude in enumerate(amplitudes):
                    if stim_type in [0,'0']:
                        exp = 0
                        stim_type = 'g'
                        electrode_temp = 0
                    else:
                        electrode_temp = electrode
                        exp = 3
                    if len(networks) == 3:
                        n = i+j+k+l
                        im = plot_image(axs[n], images[n,:,:], **get_params(exp, electrode_temp, stim_type, amplitude, networks))
                    elif len(networks) == 2:
                        for m,network in enumerate(networks):
                            im = plot_image(axs[m], images[m,:,:], **get_params(exp, electrode_temp, stim_type, amplitude, [network]))

    for i,label in enumerate(labels):
        axs[i].set_title(label, y=1.1, size=30)  
    handles, labels = axs[0].get_legend_handles_labels()
    fig.tight_layout(h_pad=0.0,w_pad=.5)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=18, shrink=0.60)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, 0), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()


def plot_1x2(exps, electrodes, stim_types, amplitudes, networks, save=None, labels=None, radius=200):

    ### Initialisation
    images = np.zeros((2,radius*2,radius*2))
    fig, axs = plt.subplots(1,2, sharex=True, sharey=True, figsize=(11,7))

    for i,exp in enumerate(exps):
        for j,electrode in enumerate(electrodes):
            for k,stim_type in enumerate(stim_types):
                for l,amplitude in enumerate(amplitudes):
                    # if stim_type in [0,'0']:
                    #     exp = 0
                    #     stim_type = 'g'
                    #     electrode_temp = 0
                    # else:
                    #     electrode_temp = electrode
                    #     exp = 3
                    # if len(networks) == 3:
                    n = i+j+k+l
                    images[n,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
                    # elif len(networks) == 2:
                    #     for m,network in enumerate(networks):
                    #         images[m,:,:] = get_image(**get_params(exp, electrode_temp, stim_type, amplitude, [network]))

    images = np.array(images)/np.nanmax(images)

    for i,exp in enumerate(exps):
        for j,electrode in enumerate(electrodes):
            for k,stim_type in enumerate(stim_types):
                for l,amplitude in enumerate(amplitudes):
                    # if stim_type in [0,'0']:
                    #     exp = 0
                    #     stim_type = 'g'
                    #     electrode_temp = 0
                    # else:
                    #     electrode_temp = electrode
                    #     exp = 3
                    # if len(networks) == 3:
                    n = i+j+k+l
                    im = plot_image(axs[n], images[n,:,:], **get_params(exp, electrode, stim_type, amplitude, networks))
                    # elif len(networks) == 2:
                    #     for m,network in enumerate(networks):
                    #         im = plot_image(axs[m], images[m,:,:], **get_params(exp, electrode_temp, stim_type, amplitude, [network]))

    for i,label in enumerate(labels):
        axs[i].set_title(label, y=1.1, size=30)  
    handles, labels = axs[0].get_legend_handles_labels()
    fig.tight_layout(h_pad=0.0,w_pad=.5)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=18, shrink=0.6)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, 0), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()


def plot_depth(exp, electrode, stim_type, amplitude, networks, save=None, labels=None, depths=None):

    ### Initialisation
    radius = 400 if exp in ['1',1] else 200
    images = np.zeros((len(depths),radius*2,radius*2))
    fig, axs = plt.subplots(1,len(depths), sharex=True, sharey=True, figsize=(16,7))


    images[0,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']))
    for i in range(1,len(depths)):
        depth = [np.min(depths[i-1:i+1]), np.max(depths[i-1:i+1])]
        images[i,:,:] = get_image(**get_params(exp, electrode, stim_type, amplitude, ['0','1','2']), depth=depth)

    images = np.array(images)/np.nanmax(images)

    im = plot_image(axs[0], images[0,:,:], **get_params(exp, electrode, stim_type, amplitude, networks))
    for i in range(1,len(depths)):
        depth = [np.min(depths[i-1:i+1]), np.max(depths[i-1:i+1])]
        im = plot_image(axs[i], images[i,:,:], **get_params(exp, electrode, stim_type, amplitude, networks), depth=depth)
            

    if labels is not None:           
        for i,label in enumerate(labels):
            axs[i].set_title(label, y=1.1, size=30)  
    handles, labels = axs[0].get_legend_handles_labels()
    fig.tight_layout(h_pad=0.0,w_pad=.5)
    fig.colorbar(im, ax=axs, orientation='vertical', aspect=18, shrink=0.6)
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.45, 0), framealpha=0)
    if save is not None:
        fig.savefig(save+'.png', dpi=400, bbox_inches='tight')
        fig.savefig(save+'.svg', bbox_inches='tight')
    else:
        plt.show()
