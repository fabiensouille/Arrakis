import os
import numpy as np

def load_res(resname="RES", record=0):
    """ 
    Load res file 
    """
    # Define result file name depending on record
    if record == -1:
        resfile = os.path.join('RESU', resname + 'fin.dat')
    elif record == -2:
        resfile = resname
    else:
        numb = str("{}".format(record))
        numb = "0"*(4-len(numb)) + numb
        resfile = os.path.join('RESU', resname + numb + '.dat')
        
    # read global parameters
    f = open(resfile, 'r')
    line = f.readline()
    line = line.split(',')
    time = float(line[0]) 
    nx = int(line[1])

    # read header
    line = f.readline()
    variables = line.split(',')
    variables = [s.strip() for s in variables if (s!='0')and(s!='0\n')]
    f.close()

    # read tables
    res = np.loadtxt(resfile, delimiter=',', skiprows=2)
    
    return time, res, variables
