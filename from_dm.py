#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:24:39 2020

@author: nico
"""
from pyscf import gto
from pyscf.dft.numint import eval_ao, eval_rho
import dmtools.dmtools as dmt
import cubetools as ct
from geomtools.geom import geom
import numpy as np

basfile="sto-6g.1.nwchem"
geomfile="acrolein.xyz"   
dmfile="frag0_HF_dens_ME.txt"
Np_vect=np.array([75,75,75])
Vect_M=np.array([[0.2,0,0],[0,0.2,0],[0,0,0.2]])
origin=np.array([-11.0,-6.0,-4.0])
geometry = geom.from_xyz(geomfile)
geomstring =geometry.__str__()
with open(basfile, "r") as f:
    ibasis = f.read()
mol = gto.M(atom = geomstring, basis = ibasis)
dm_obj = dmt.DM.from_dmfile(dmfile)
dm = dm_obj.get_dm_full()
values = np.zeros(list(Np_vect))
geometry.change_coord_unit("au")
cube = ct.cubefile(geometry,Np_vect,Vect_M,values,origin,comment="test\ntest\n")
points=cube.get_coordlist()
ao_mol = eval_ao(mol, points, deriv=0)
cubevals= eval_rho(mol, ao_mol, dm, xctype='LDA')
cube.values = cubevals.reshape(*Np_vect)
cube.write("cube_from_dm.cube")


