[Back to ToC](/docs/manual/README.md)

# Calculating extracellular potentials with COMSOL

In order to calculate the extracellular potentials inside the tissue that arise from electrical stimulation, we require a FEM model of tissue on which the right boundary conditions are imposed. The whole process of creating the model geometry, generating a mesh, assigning materials, choosing physics, imposing boundary conditions... can be done in COMSOL.

## Study

We are interested in the spatial and temporal behaviour of the extracellular potentials in response to certain current injections. We could represent these space- and time-dependent potentials in a matrix. There are $N$ rows representing the FEM mesh nodes, and $T$ columns representing the timepoints. Let's call this matrix $\bf{S}$ for solution.

$$ \bf{S} = 
\begin{bmatrix}
V_{X_1,t_1} & V_{X_1,t_2} & \cdots & V_{X_1,t_T} \\
V_{X_2,t_1} & V_{X_2,t_2} & \cdots & V_{X_2,t_T} \\
\vdots      & \vdots      & \ddots & \vdots          \\
V_{X_N,t_1} & V_{X_N,t_2} & \cdots & V_{X_N,t_T} \\
\end{bmatrix}
$$

### 1. One time-dependent study

In the most general case, we need to solve the FEM at every timepoint (time-dependent study), meaning the full solution is described by $ N \times T $ values. Using the output of such a study is supported in the comsol module, but is not very efficient (w.r.t computation time and storage) and should thus only be used if the method in the paragraph below is not an option. 


### 2. One stationary study

Thanks to the quasi-static approximation, the FEM solution is linear w.r.t. the injected current(s). In cases that are not too complex, i.e. where the same current profile (but with possibly different amplitudes) is used for all electrodes, the FEM solutions at different time points are a function of just the current amplitude, meaning they are linearly dependent. Such a matrix S is of rank 1 and can be written as the outer product of two vectors.

$$ \bf{S} = 
\begin{bmatrix}
V_{X_1,t_1} & V_{X_1,t_2} & \cdots & V_{X_1,t_T} \\
V_{X_2,t_1} & V_{X_2,t_2} & \cdots & V_{X_2,t_T} \\
\vdots      & \vdots      & \ddots & \vdots          \\
V_{X_N,t_1} & V_{X_N,t_2} & \cdots & V_{X_N,t_T} \\
\end{bmatrix}
= \begin{bmatrix}
V_{X_1}A_{t_1} & V_{X_1}A_{t_2} & \cdots & V_{X_1}A_{t_T} \\
V_{X_2}A_{t_1} & V_{X_2}A_{t_2} & \cdots & V_{X_2}A_{t_T} \\
\vdots      & \vdots      & \ddots & \vdots          \\
V_{X_N}A_{t_1} & V_{X_N}A_{t_2} & \cdots & V_{X_N}A_{t_T} \\
\end{bmatrix} 
= \vec{V}_X \otimes \vec{A}_t
$$

Here, the full solution can be described by the FEM solution at one time point ($\vec{V_X}$) and a time-dependent scaling factor (i.e. the current profile $\vec{A_t}$), resulting in only $ N + T $ values. This is a supported (and recommended) way to define the extracellular potentials in the comsol module of BioNet.


### 3. Multiple stationary studies

Because of the linearity of the solutions, the full solution $\bf{S}$ can also be defined as the superposition (i.e. linear combination) of the solutions $\bf{S}_i$ where each electrode is active by itself.

$$ \bf{S} = \sum_i \bf{S}_i = \sum_i \vec{V}\_{X,i} {\otimes} \vec{A}\_{t,i} $$

When only one electrode is active, the solution can always be decomposed into a spatial component and a temporal component as in the paragraph above. Doing this decomposition for each electrode separately and adding the solutions, only requires the FEM to be solved once for each electrode, and requires $ N_{electrodes} \times (N + T) $ values to store the full solution. In almost every case, this would be faster and more storage-efficient than the method described in the first paragraph. However, this is not (yet) supported in the comsol module.


## Output

After a solution has been calculated, it can be exported with Results>Export>Data.

- File type: Text
- Points to evaluate in: Take from dataset
- Data format: Spreadsheet

This will generate a .txt file with a bunch of header rows (starting with %), and then at least 4 space-separated columns. The first three columns are the x-, y-, and -coordinate, where every row defines the 3D-coordinates of one of the mesh nodes.

Depending on whether simulation was stationary or time-dependent, there will be either one or multiple extra columns.
- Stationary: The 4th column described the potential at each point. This column is essentially $ \vec{V_X} $.
- Time-dependent: Every columns from the 4th on contains the voltage profile at one timepoint T, like a vector $\begin{bmatrix} V_{X_0,t_T} & V_{X_1,t_T} & \cdots & V_{X_N,t_T}\end{bmatrix}^T$ that corresponds to a column of matrix $S$.

Configuring the comsol input for BMTK in the config.json file, will look something like this

```json
    "Extracellular_Stim": {
        "input_type": "lfp",
        "node_set": "all",
        "module": "comsol",
        "comsol_file": "$STIM_DIR/exp3/02-.txt",
        "waveform": "$STIM_DIR/waveform.csv",
        "amplitude": 10,
        "ip_method": "L"
    }
```
- input_type - Always "lfp".
- node_set - Used to filter which cells receive the input, but here it probably does not make sense to use anything besides "all".
- module - Always "comsol".
- comsol_file - /path/to/comsol.txt from stationary or time-dependent study
- waveform - /path/to/waveform.csv. Only supply this with stationary comsol.txt, remove the entire line for time-dependent comsol.txt.
- amplitude - Multiplication factor for waveform. If the amplitudes in waveform.csv are normalised to [0,1], this can be used to set the current amplitude. Defaults to 1. 
- ip_method - "NN" for nearest neighbour, "L" for linear interpolation (only for stationary comsol.txt). Defaults to "NN".
