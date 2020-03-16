#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:24:39 2020

@author: nico
"""
import cubetools.cubetools as ct
from geomtools.geom import geom
import numpy as np

basfile="sto-6g.1.nwchem"
geomfile="acrolein.xyz"   
dmfile="frag0_HF_dens_ME.txt"
Np_vect=np.array([75,75,75])
Vect_M=np.array([[0.2,0,0],[0,0.2,0],[0,0,0.2]])
origin=np.array([-11.0,-6.0,-4.0])
g = geom.from_xyz(geomfile)

cf_empty= ct.read_cubefile("N_HF.cube",empty_values=True) #how to get an empty cubefile obj from another cubefile
cf_empty2 = ct.cubefile.from_grid_specs(g=g,Np_vect=Np_vect, Vect_M=Vect_M, o=origin) #how to get an empty cubefile obj by defining the grid

cf_empty2.get_cubevals_from_dm(dmfile=dmfile,basfile=basfile) #replacing the empty values with those obtained from DM
fromdm = ct.get_cube_from_dm(dmfile=dmfile,basfile=basfile,g=g,Np_vect=Np_vect,Vect_M=Vect_M,o=origin) # same as the two previous step but in one line



