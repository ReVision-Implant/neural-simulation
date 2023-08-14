import os
import pandas as pd
import numpy as np

dir_out = './examples/bio_components/stimulations/exp0'
dir_in = './examples/bio_components/stimulations/exp0'
for f in os.listdir(dir_in):
    # if f[1] not in ['1', '2', '3']:
    #     continue
    path = os.path.join(dir_in,f)
    header = pd.read_csv(path, sep=';', nrows=9,header=None)
    body = pd.read_csv(path, sep='\s+',skiprows=9, header=None, engine='python')

    header[1] = [' ',' ',' ',' ',' ',' ',' ',' ',' ']
    header[2] = [' ',' ',' ',' ',' ',' ',' ',' ',' ']
    header[3] = [' ',' ',' ',' ',' ',' ',' ',' ',' ']

    body[0] -= 75
    # body[1] = 800-body[1]
    body[2] -= 91

    df = pd.concat([header, body], ignore_index=True)
 
    if os.path.exists(os.path.join(dir_out,f)):
        os.remove(os.path.join(dir_out,f))
    with open(os.path.join(dir_out,f), 'a') as file:
        df_string = df.to_string(index=False, header=False, justify='left')
        file.write(df_string)

    header = pd.read_csv(os.path.join(dir_out,f), sep="\s{3,}", header=None, skiprows=8, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
    header[0] = header[0][2:]                               # remove '% ' before first column name
    if header[3][3] == 'V':
        unit = 1000
    elif header[3][3] == 'm':
        unit = 1
    for i,col in enumerate(header):                         # remove superfluous characters before actual time value
        if col[0] == "V":
            header[i] = 0
    timepoints = np.array(header[3:], dtype=float)    # create array of timepoints  

    # load data in COMSOL output .txt file.  
    comsol = pd.read_csv(os.path.join(dir_out,f), sep="\s+", header=None, skiprows=9, names=header)

    print(comsol)