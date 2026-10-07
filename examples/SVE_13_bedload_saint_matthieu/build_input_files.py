import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

if __name__ == '__main__':
       
    # run init 
    print("~~> run init")
    os.system("arrakis.py saint_matthieu_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
