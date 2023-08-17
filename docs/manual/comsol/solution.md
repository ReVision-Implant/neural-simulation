[Back to ToC](/docs/manual/README.md)

# Calculating extracellular potentials with COMSOL

In order to calculate the extracellular potentials inside the tissue that arise from electrical stimulation, we require a FEM model of tissue on which the right boundary conditions are imposed. The whole process of creating the model geometry, generating a mesh, assigning materials, choosing physics, imposing boundary conditions... can be done in COMSOL.

## Study

We are interested in the spatial and temporal behaviour of the extracellular potentials in response to certain current injections. We could represent these space- and time-dependent potentials in a 2D matrix. One dimension represents the spatial component and holds all domain nodes, and the other dimensions represents the temporal component and holds all time points. Let's call this matrix S for solution.

$$ S = 
\begin{bmatrix}
V_{X_0,t_0} & V_{X_0,t_1} & \cdots & V_{X_1,t_{max}} \\
V_{X_1,t_0} & V_{X_1,t_1} & \cdots & V_{X_1,t_{max}} \\
\vdots      & \vdots      & \ddots & \vdots          \\
V_{X_N,t_0} & V_{X_N,t_1} & \cdots & V_{X_N,t_{max}} \\
\end{bmatrix}
$$

Thanks to the quasi-static approximation, the FEM solutions are linear. In stimulation cases that are not too complex, i.e. where only a single current profile (amplitudes can be different) is used for all electrodes, the voltage profile in the tissue will remain the same and the voltage at every point will only vary proportionally to the current amplitude. In such case, the matrix can be written as the outer product of two vectors ($\vec{V_X}$).

$$ S = 
\begin{bmatrix}
V_{X_0}A_{t_0} & V_{X_0}A_{t_1} & \cdots & V_{X_0}A_{t_{max}} \\
V_{X_1}A_{t_0} & V_{X_1}A_{t_1} & \cdots & V_{X_1}A_{t_{max}} \\
\vdots      & \vdots      & \ddots & \vdots          \\
V_{X_N}A_{t_0} & V_{X_N}A_{t_1} & \cdots & V_{X_N}A_{t_{max}} \\
\end{bmatrix} 
= \vec{V_X} \otimes \vec{A_t}
$$

As such, we only need to solve the FEM at one time point (stationary study) and we only need $ N + t_{max} $ values to describe the time-dependent solution of every point inside the tissue. This is a supported way to define the extracellular potentials in the comsol module of BioNet.

For very complex stimulation algorithms, where multiple current profiles are used (e.g. using slightly different current pulses timings on different electrodes), you will need to solve the FEM at every timepoint (time-dependent study), and we need $ N \times t_{max} $ values to describe the time-dependent solution. This is also supported in the comsol module . 

In theory, this last situation could also be described by solving the FEM model at one time point for each electrode (one stationary study per electrode), and then using the superposition principle to combine the results as if all electrodes were on simultaneously. This requires the FEM to be solved once for each electrode, and requires $ N_{electrodes} \times (N + t_{max}) $ values to store the full solution. This would be faster and more storage-efficient than the method described in the previous paragraph. However, this is not (yet) supported in the comsol module.


## Output

After a solution has been calculated, it can be exported with Results>Export>Data.

- File type: Text
- Points to evaluate in: Take from dataset
- Data format: Spreadsheet

This will generate a .txt file with a bunch of header rows, and then at least 4 space-separated columns. The first three columns are always the X-, Y-, and Z-coordinate, where every row defines the 3D-coordinates of one of the mesh nodes.

Depending on whether simulation was stationary or time-dependent, there will be either one or multiple extra columns.
- Stationary: The 4th column described the potential at each point. This column is essentially $ \vec{V_X} $ in the equations above.
- Time-dependent: Every columns from the 4th on contains the voltage profile at one timepoint T, like a vector $\begin{bmatrix} V_{X_0,t_T} & V_{X_1,t_T} & \cdots & V_{X_N,t_T}\end{bmatrix}^T$ that is a column of matrix $S$ above.