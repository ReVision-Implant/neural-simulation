import os
import pandas as pd
import numpy as np

def translate_comsol(dir_in, dir_out, X, Y, Z):
    """This function translates the node coordinates of a comsol.txt file.
    Can be used to reuse comsol solutions with different coordinates of the BMTK network.
    The function iterates over all files in dir_in, 
    translates them according to the lambda expressions in x,y,z, 
    and saves them in dir_out.

    :param dir_in: path/to/dir of original comsol.txt files
    :type dir_in: str
    :param dir_out: path/to/dir of translated comsol.txt files
    :type dir_out: str
    :param x: x_tr = X(x)
    :type x: lambda
    :param y: y_tr = Y(y)
    :type y: lambda
    :param z: z_tr = Z(z)
    :type z: lambda
    """    

    for f in os.listdir(dir_in):
        if 
        path = os.path.join(dir_in,f)
        header = pd.read_csv(path, sep=';', nrows=9, header=None)
        body = pd.read_csv(path, sep='\s+',skiprows=9, header=None, engine='python')

        # Add empty columns 
        header[1] = [' ',' ',' ',' ',' ',' ',' ',' ',' ']
        header[2] = [' ',' ',' ',' ',' ',' ',' ',' ',' ']
        header[3] = [' ',' ',' ',' ',' ',' ',' ',' ',' ']

        body[0] = X(body[0])
        body[1] = Y(body[1]) 
        body[2] = Z(body[2])

        df = pd.concat([header, body], ignore_index=True)
    
        if os.path.exists(os.path.join(dir_out,f)):
            os.remove(os.path.join(dir_out,f))
        with open(os.path.join(dir_out,f), 'a') as file:
            df_string = header.to_string(index=False, header=False, justify='left')
            df_string = df_string.lstrip()
            file.write(df_string)

        # header = pd.read_csv(os.path.join(dir_out,f), sep="\s{3,}", header=None, skiprows=8, nrows=1, engine='python').to_numpy()[0] # load header row of .txt file
        # header[0] = header[0][2:]                               # remove '% ' before first column name
        # if header[3][3] == 'V':
        #     unit = 1000
        # elif header[3][3] == 'm':
        #     unit = 1
        # for i,col in enumerate(header):                         # remove superfluous characters before actual time value
        #     if col[0] == "V":
        #         header[i] = 0
        # timepoints = np.array(header[3:], dtype=float)    # create array of timepoints  

        # # load data in COMSOL output .txt file.  
        # comsol = pd.read_csv(os.path.join(dir_out,f), sep="\s+", header=None, skiprows=9, names=header)

        # print(comsol)

translate_comsol(
    dir_in = './components/stimulations/exp0',
    dir_out = './components/stimulations/exp0_test',
    X = lambda x:x+100,
    Y = lambda y:y-100,
    Z = lambda z:z
)
