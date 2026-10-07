import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def zb(x):
    if x >= 4. and x < 280.:
      if x >= 12. and x < 212.:
        gx = 0.02*np.sin(0.04*np.pi*(x-12.))
      else:
        gx = 0.
      return (0.2/276.)*(x-4.) + gx
    elif x >= 280. and x < 556.:
      return 0.2 - (0.2/276.)*(x-280.)
    else:
      return 0.

def h(x):
    return max(0.7 - zb(x), 0.)

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=140,
        xa=0.,
        xb=560.,
        geo_file="MESH.geo",
        zb_funct=zb,
        verbose=False)

    # build initial condition for init:
    print("~~> create initial condition for init")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h, 
        u=1.,
        verbose=False)
        
    # run init 
    print("~~> run init")
    os.system("arrakis.py travelling_jump_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
