#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:29:10 2019

@author: nico
"""
from geomtools import geom
import numpy as np
import re,sys

class cubefile:
    def __init__(self, geom, Np_vect, Vect_M, values, origin=np.zeros(3), comment="\n \n"):
#    def __init__(self, geom, Np_vect, Vect_M, values, origin=np.zeros(3)):
        self.geom=geom
        self.Np_vect=Np_vect
        self.Vect_M=Vect_M
        self.values=values
        self.origin=origin
        self.comment=comment
    
    def from_file(fname):
        return read_cubefile(fname)
    def write(self, fname, digits=5):
        import mendeleev as md
        with open(fname, "w") as out:
            out.write(self.comment)
            out.write("{:5}{:12.6f}{:12.6f}{:12.6f}\n".format(len(self.geom.atoms),self.origin[0],self.origin[1],self.origin[2]))
            for i in range(3):
                out.write("{:5}{:12.6f}{:12.6f}{:12.6f}\n".format(self.Np_vect[i],self.Vect_M[i][0],self.Vect_M[i][1],self.Vect_M[i][2]))
            for i in range(len(self.geom.atoms)):
#                out.write(self.geom.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}'.format(g.coords(which)[i][0],g.coords(which)[i][1],g.coords(which)[i][2],w=spacing+decs+1, p=decs)+"\n")
                out.write('{:5}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n'.format(md.element(self.geom.atoms[i]).atomic_number,md.element(self.geom.atoms[i]).atomic_number,self.geom.inp_coords[i][0],self.geom.inp_coords[i][1],self.geom.inp_coords[i][2]))
            for x in range(self.Np_vect[0]):
                for y in range(self.Np_vect[1]):
                    for z in range(self.Np_vect[2]):
                        if (z+1) % 6 == 0 or z == self.Np_vect[2]-1:
                            out.write("{:{w}.{dd}e}\n".format(self.values[x][y][z],w=8+digits,dd=digits))
                        else:
                            out.write("{:{w}.{dd}e}".format(self.values[x][y][z],w=8+digits,dd=digits))
    def add(self, cf):
        return add(self, cf)
    def subtract(self, cf):
        return subtract(self, cf)
    def times(self, factor):
        self.values=np.multiply(factor, self.values)
    def trim(self, last=True):
        if last==True:
            self.values=self.values[:-1,:-1,:-1]
        else:
            self.values=self.values[1:,1:,1:]
        self.Np_vect-=1
    def RSR(self,cf):
        return RSR(self,cf)
    def num(self,cf):
        return num(self,cf)
    def den(self,cf):
        return den(self,cf)
    def get_coords(self,arr_pos): #nb. single point should be still given as array [[x],[y],[z]]
        all_coords=np.zeros(arr_pos.shape)
        for n,i in enumerate(arr_pos):
            coords=np.copy(self.origin)
            for j in range(3):
                coords+=i[j]*self.Vect_M[j]
            all_coords[n]=coords
        return all_coords
    def get_pointlist(self):
        return np.transpose(np.where(self.values!=np.nan)) #trick to get it always verified
    def get_coordlist(self):
        return self.get_coords(self.get_pointlist())
    def min_dist(self,point):
        d=np.zeros(len(self.geom.atoms))
        for n,i in enumerate(self.geom.inp_coords):
            d[n]=np.linalg.norm(i-point)
        return d.min()
        
def read_cubefile(fname):
#fname="test.cube"
    import mendeleev as md
    pnums = []
    V_m=[]
    cmmnt=""
    natm = 1
    positions = []
    atms=[]
    with open(fname) as f:
        for j,line in enumerate(f):
            if j < 2:
                cmmnt+=line
            if j == 2:
                splt =line.split()
                natm = int(splt[0])
                o = np.asarray(list(map(float,splt[1:4])))
            if j > 2 and j < 6:
                splt =line.split()
                pnums.append(int(splt[0]))
                V_m.append(np.asarray(list(map(float,splt[1:4]))))
            if j >= 6 and j < natm+6:
                splt=line.split()
                atms.append(md.element(int(splt[0])).symbol)
                positions.append(np.asarray(list(map(float,splt[2:5]))))
            if j == natm+6:
                break
    Np_vect = np.asarray(list(map(int,pnums)))
    Vect_M=np.asarray(V_m)
    coords=np.asarray(positions)
    atoms=np.asarray(atms)
    g=geom.geom(atoms,coords, coord_unit="au")
    with open(fname) as f:
        cube = f.read()
    cubelist=[]
    pattern2 = '(-?\d+\.?\d+[Ee]\ *[-+]?\ *\d+\s*){' + re.escape(str(pnums[2])) + '}'
    cubelist.append([])
    xcnt = 0
    ycnt = 0
    for xy,m in enumerate(re.finditer(pattern2,cube)):
        if xy % Np_vect[1] == 0 and xy > 0:
            cubelist[xcnt]=np.asarray(cubelist[xcnt])
            xcnt += 1 if xy>0 else 0
            ycnt = 0
            cubelist.append([]) # append empty list for x
        s = m.group().split()
        s = np.asarray(list(map(float,s)))
        cubelist[xcnt].append(s)
        ycnt += 1
    cubelist[-1]=np.asarray(cubelist[-1])
    cubearray=np.asarray(cubelist)
    return cubefile(g,Np_vect,Vect_M,cubearray,origin=o,comment=cmmnt)

def check_same_grid(cf1,cf2,orig_thresh=0.00001):
    if not (cf1.Np_vect==cf2.Np_vect).all():
        print("Not the same number of points per direction!!")
        return False
    elif not (cf1.Vect_M==cf2.Vect_M).all():
        print("Not the same grid vectors!!")
        return False
    elif not ((cf1.origin-cf2.origin) < orig_thresh).all():
        print("Not the same origin!!")
        return False
    elif cf1.values.shape!=cf2.values.shape:
        print("Not the same number of total points!")
        return False
    else:
        return True

def subtract(cf1,cf2,comment="subtracted with cubetools\n here it goes\n",geom="",Continue=False):
    if comment !="subtracted with cubetools\n here it goes\n":
        if len(comment.split("\n"))==3:
            pass
        elif len(comment.split("\n")) < 3:
            comment+="\n"*(3-len(comment.split("\n")))
        elif len(comment.split("\n")) > 3:
            print("Please retry with a two- or fewer-lines comment")
            sys.exit()
    if geom=="":
        geom=cf1.geom
    if check_same_grid(cf1,cf2):
        Np_vect=cf1.Np_vect
        Vect_M=cf1.Vect_M
        origin=cf1.origin
    else:
        if Continue:
            print("Not the same grid! creating mock-object and continuing")
            Np_vect=cf1.Np_vect
            Vect_M=cf1.Vect_M
            origin=cf1.origin
            return cubefile(geom,Np_vect,Vect_M,np.zeros(cf1.values.shape),origin,comment=comment)
        else:
            print("Not the same grid! Exiting!")
            sys.exit()
    return cubefile(geom,Np_vect,Vect_M,np.subtract(cf1.values,cf2.values),origin,comment=comment)

def add(cf1,cf2,comment="subtracted with cubetools\n here it goes\n",geom="",Continue=False):
    if comment !="subtracted with cubetools\n here it goes\n":
        if len(comment.split("\n"))==3:
            pass
        elif len(comment.split("\n")) < 3:
            comment+="\n"*(3-len(comment.split("\n")))
        elif len(comment.split("\n")) > 3:
            print("Please retry with a two- or fewer-lines comment")
            sys.exit()
    if geom=="":
        geom=cf1.geom
    if check_same_grid(cf1,cf2):
        Np_vect=cf1.Np_vect
        Vect_M=cf1.Vect_M
        origin=cf1.origin

        if Continue:
            print("Not the same grid! creating mock-object and continuing")
            Np_vect=cf1.Np_vect
            Vect_M=cf1.Vect_M
            origin=cf1.origin
            return cubefile(geom,Np_vect,Vect_M,np.zeros(cf1.values.shape),origin,comment=comment)
        else:
            print("Not the same grid! Exiting!")
            sys.exit()
    return cubefile(geom,Np_vect,Vect_M,np.add(cf1.values,cf2.values),origin,comment=comment)

def RSR(cf1,cf2,Continue=False):
    if not check_same_grid(cf1,cf2):
        if Continue:
            print("Not the same grid, but continuing!Returning 2!")
            return 2.0
        else:
            print("Not the same grid! Exiting!")
            sys.exit()
    else:            
        num=(abs(np.subtract(cf1.values,cf2.values))).sum()
        den=(abs(np.add(cf1.values,cf2.values))).sum()
    return np.divide(num,den)

def num(cf1,cf2,Continue=False):
    if not check_same_grid(cf1,cf2):
        if Continue:
            print("Not the same grid, but continuing!Returning 2!")
            return 2.0
        else:
            print("Not the same grid! Exiting!")
            sys.exit()
    else:            
        num=(abs(np.subtract(cf1.values,cf2.values))).sum()
    return num

def den(cf1,cf2,Continue=False):
    if not check_same_grid(cf1,cf2):
        if Continue:
            print("Not the same grid, but continuing!Returning 2!")
            return 2.0
        else:
            print("Not the same grid! Exiting!")
            sys.exit()
    else:            
        den=(abs(np.add(cf1.values,cf2.values))).sum()
    return den
#def RSR(cf1,cf2,Continue=False):
#    if not check_same_grid(cf1,cf2):
#        if Continue:
#            print("Not the same grid, but continuing!Returning 2!")
#            return 2.0
#        else:
#            print("Not the same grid! Exiting!")
#            sys.exit()
#    else:
#        diff=np.subtract(cf1.values,cf2.values)
#        adiff=abs(diff)            
#        num=adiff.sum()
#        add=np.add(cf1.values,cf2.values)
#        aadd=abs(add)
#        if not (add==aadd).all():
#            sys.exit()
#        den=aadd.sum()
#        print("MinMax diff", diff.min(), diff.max())
##        print("MinMax adiff", adiff.min(), adiff.max())
#        print("MinMax add", add.min(), add.max())
#    return num,den
    


    