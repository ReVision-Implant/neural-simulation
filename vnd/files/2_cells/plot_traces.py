import numpy as np
import h5py
import matplotlib.pyplot as plt


def plot_results(node_ids, pop_name, labels=None):
    with h5py.File('output/v_report.h5', 'r') as h5:
        fig, axes = plt.subplots(len(node_ids))
        for i, node_id in enumerate(node_ids):
            report_grp = h5['/report'][pop_name]
            index_pointer_arr = report_grp['mapping/index_pointer'][()]
            node_idx = int(report_grp['mapping/node_ids'][node_id][()])
            idx_beg, idx_end = index_pointer_arr[node_idx], index_pointer_arr[node_idx+1]
            seg_types_arr = report_grp['mapping/seg_types'][idx_beg:idx_end]
            vm_data = report_grp['data'][:, idx_beg:idx_end]

            ts_arr = report_grp['mapping/time'][()]
            ts_beg, ts_end, ts_dt = ts_arr[0], ts_arr[1], ts_arr[2]
            ts = np.linspace(ts_beg, ts_end, num=int((ts_end-ts_beg)/ts_dt))

            # plot mean soma trace
            soma_idxs = np.argwhere(seg_types_arr == 1).flatten()
            mean_soma_vm = np.mean(vm_data[:, soma_idxs], axis=1)
            axes[i].plot(ts, mean_soma_vm, label='soma')

            # plot mean basal dend trace
            dend_idxs = np.argwhere(seg_types_arr == 3).flatten()
            mean_dend_vm = np.mean(vm_data[:, dend_idxs], axis=1)
            axes[i].plot(ts, mean_dend_vm, label='dend')

            # plot mean apical dend traces
            apic_idxs = np.argwhere(seg_types_arr == 4).flatten()
            if len(apic_idxs) > 0:
                apic_dend_vm = np.mean(vm_data[:, apic_idxs], axis=1)
                axes[i].plot(ts, apic_dend_vm, label='apic')

            axes[i].set_xlabel('time (ms)')
            axes[i].set_ylabel(report_grp['data'].attrs['variable'])
            axes[i].legend()
            if i + 1 < len(node_ids):
                axes[i].get_xaxis().set_visible(False)

            if labels is not None:
                axes[i].title.set_text(labels[i])
            else:
                axes[i].title.set_text(f'node_id {node_id}')

        plt.savefig('v_report_seg_avgs.png')
        plt.show()


if __name__ == '__main__':
    plot_results(node_ids=[0, 1], pop_name='v1', labels=['Scnn1a', 'Pvalb'])
