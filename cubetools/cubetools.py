#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:29:10 2019

@author: nico
"""
from geomtools.geom import geom
import numpy as np
import re

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class gridError(Error):
    """Mismatch in the grids.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message = "Mismatch in grids!!"):
        self.message = message
        print(self.message)
        
class otherError(Error):
    """Mismatch in the units.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        print(self.message)
        
class cubefile:
    """
    Note
    ----
    main class in this molecule. Contains all information of a cubefile.
    """
    def __init__(self, geom, Np_vect, Vect_M, values, origin=np.zeros(3), comment="\n \n"):
        """
        geom: geom
            geometry object(geomtools) 
        Np_vect: array(3,)
            number of points along x,y,z
        Vect_M: array(3,3)
            grid vectors. Generally diagonal
        values: array(Np_vect[0],Np_vect[1],Np_vect[2])
            values of the function on the grid
        origin: array(3,)
            origin of the grid. default is 0,0,0
        comment: str
            optional comment.
        """
        self.geom = geom
        self.Np_vect = Np_vect
        self.Vect_M = Vect_M
        self.values = values
        self.origin = origin
        self.comment = comment
    
    def copy(self):
        """
        Returns
        -------
        cubefile
            a copy of self
        """
        import copy as c
        return c.deepcopy(self)
    
    def from_file(fname, emptyval=False):
        """
        Note
        ----
        Allows quick reading of a cubefile.
        
        Parameters
        ----------
        fname: str
            filename or path
        """
        return read_cubefile(fname)
    
    def from_grid_specs(g="", Np_vect=None, Vect_M=None,  o=None, comment="Empty function values\nOnly grid\n"):
        """
        Note
        ----
        Returns an empty-values cubefile object, useful to then fill this with values from a density matrix, or to store grid specs
        
        Parameters
        ----------
        g: geom or str
            either geom object or string of the xyz filename
        Np_vect: array(3,)
            number of points along x,y,z
        Vect_M: array(3,3)
            grid vectors. Generally diagonal
        o: array(3,)
            origin of the grid. default is 0,0,0
        comment: str
            optional comment.
        
        Returns
        -------
        cubefile obj
            Cubefile with desired grid but empty value array
        """
        return emptyval_cube_from_specs(g=g, Np_vect=Np_vect, Vect_M=Vect_M,  o=o, comment=comment)
    
    def get_cubevals_from_dm(self, dmfile="", basfile="", g="", comment="Obtained from density matrix in file: \n {dmfile} \n"):
        """
        Note
        ----
        Allows to obtain values from a density matrix.
        
        Parameters
        ----------
        dmfile: str or list
            density matrix file, list or "file1,file2"
        basfile: str or list
            basis .nwchem file, list or "file1.nwchem,file2.nwchem"
        g: str or geom obj, or [str,str] or [obj,obj]
            geometry or geometry file which overrides self.geom. Compulsory only for dual basis. 
            For dual basis: "geomA.xyz, geomB.xyz",  ["geomA.xyz"," geomB.xyz"], or [geomA, geomB]
        comment: str
            default writes dmfile, use empty string to leave self.cube
        Sets
        ----
        self.values
            obtained from density matrix
        self.comment
            user-chosen, default is corresponding density matrix
        """
        basfile = list(basfile) if type(basfile) == tuple else basfile  # no tuples because immutable
        g = list(g) if type(g) == tuple else g  # no tuples because immutable
        from pyscf import gto
        from pyscf.dft.numint import eval_ao, eval_rho
        import dmtools.dmtools as dmt
        if dmfile == "":
            dmfile = input("You must specify the density matrix file. Can have header or not, can be alpha and beta or not. Please type the filename\n") #TODO check if it reads a single DM as 2*a or a.
        if basfile == "":
            import glob as gl
            files = gl.glob("*.nwchem")
            if len(files) == 0:
                raise FileNotFoundError("no *.nwchem file! basis must be in .nwchem!")
            elif len(files) == 1:
                basfile = files[0]
                print("basfile not specified, using {} as a single basis".format(basfile))
            else: 
                basfile = input("Too many *.nwchem. For single basis, type the one to use. For dual basis, type '[basisA.nwchem],[basisB.nwchem]'\n")
        basfile = [i.strip() for i in basfile.split(",")] if "," in basfile else basfile  # turn string to list if dual basis
        if type(basfile) == str:  # then it is not dual basis
            if g != "":
                if type(g) == str:
                    g = geom.from_xyz(g)
            self.geom = g
            self.geom.change_coord_unit("Angstrom")
            
            with open(basfile, "r") as f:
                ibasis = f.read()
        if type(basfile) == list:  # then it is dual basis
            if g == "":
                g = input("From basfile I deduce dualbasis, Please type 'geomfileA,geomfileB'\n")
            if type(g) == str:
                if len(g.split(",")) == 2:
                    g = [i.strip() for i in g.split(",")]
                else:
                    raise TypeError("Cannot obtain two geometries from {}.\n From basfile I deduce dualbasis, so you must give geomA and geomB".format(g))
            for n in range(2):
                try:
                    g[n] = geom.from_xyz(g[n], identifier=n)
                except:
                    g[n].change_coord_unit("angstrom")
                    g[n].change_identifier(n)
            self.geom = g[0] + g[1]
            ###create ibasis
            for n,i in enumerate(basfile):
                with open(i,"r") as f:
                    basfile[n] = f.read()
            ibasis = {i+str(n): basfile[n] for n in range(2) for i in set(g[n].atoms)}
        geomstring = self.geom.__str__()    
        mol = gto.M(atom = geomstring, basis = ibasis)
        dm_obj = dmt.DM.from_dmfile(dmfile)
        dm = dm_obj.get_dm_full()
        self.geom.change_coord_unit("au")
        points = self.get_coordlist()
        ao_mol = eval_ao(mol, points, deriv=0)
        cubevals = eval_rho(mol, ao_mol, dm, xctype='LDA')
        self.values = cubevals.reshape(*self.Np_vect)
        if "{dmfile}" in comment:
            comment = comment.format(dmfile=dmfile)
        self.comment = comment
        
        
    def write(self, fname, digits=5):
        """
        Note
        ----
        Allows quick writing of a cubefile.
        
        Parameters
        ----------
        fname: str
            filename or path
        digits: int
            digits desired. default is 5
        """
        import mendeleev as md
        if len(self.comment.split("\n")) == 3:
            pass
        elif len(self.comment.split("\n")) < 3:
            self.comment+="\n"*(3-len(self.comment.split("\n")))
        elif len(self.comment.split("\n")) > 3:
            raise otherError("Please retry with a two- or fewer-lines comment")
        with open(fname, "w") as out:
            out.write(self.comment)
            out.write("{:5}{:12.6f}{:12.6f}{:12.6f}\n".format(len(self.geom.atoms), self.origin[0], self.origin[1], self.origin[2]))
            for i in range(3):
                out.write("{:5}{:12.6f}{:12.6f}{:12.6f}\n".format(self.Np_vect[i], self.Vect_M[i][0], self.Vect_M[i][1], self.Vect_M[i][2]))
            for i in range(len(self.geom.atoms)):
                out.write('{:5}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n'.format(md.element(self.geom.atoms[i]).atomic_number, md.element(self.geom.atoms[i]).atomic_number, self.geom.coords[i][0],
                          self.geom.coords[i][1], self.geom.coords[i][2]))
            for x in range(self.Np_vect[0]):
                for y in range(self.Np_vect[1]):
                    for z in range(self.Np_vect[2]):
                        if (z+1) % 6 == 0 or z == self.Np_vect[2]-1:
                            out.write("{:{w}.{dd}e}\n".format(self.values[x][y][z], w=8+digits, dd=digits))
                        else:
                            out.write("{:{w}.{dd}e}".format(self.values[x][y][z], w=8+digits, dd=digits))
    
    def __str__(self):
        """
        Note
        ----
        Allows easy printing
        """
        geomstring = self.geom.__str__()
        Npoints = "Npoints: x = {}, y = {}, z = {}".format(*[self.Np_vect[i] for i in range(3)])
        isdiag = (np.diag(np.diag(self.Vect_M)) == self.Vect_M).all()
        vect = "v_x = {}, vy = {}, vz = {}".format(*[np.diag(self.Vect_M)[i] for i in range(3)]) if isdiag else self.Vect_M
        orig = "origin: {}, {}, {}\n".format(*[self.origin[i] for i in range(3)])
        return "\n".join([geomstring, Npoints, vect, orig, self.comment])
    
    def __repr__(self, spacing=4, decs=6):
        """
        Note
        ----
        Allows easy calling
        """
        geomstring = self.geom.__str__()
        Npoints = "Npoints: x = {}, y = {}, z = {}".format(*[self.Np_vect[i] for i in range(3)])
        isdiag = (np.diag(np.diag(self.Vect_M)) == self.Vect_M).all()
        vect = "v_x = {}, vy = {}, vz = {}".format(*[np.diag(self.Vect_M)[i] for i in range(3)]) if isdiag else self.Vect_M
        orig = "origin: {}, {}, {}\n".format(*[self.origin[i] for i in range(3)])
        return "\n".join([geomstring, Npoints, vect, orig, self.comment])
    
    def __add__(self, other):
        """
        Note
        ----
        Allows quick addition of cubefiles.
        
        Parameters
        ----------
        other: cubefile
            cubefile to add
        
        Returns
        -------
        cubefile
            the sum
        """
        return add(self, other)
    
    def __sub__(self, other):
        """
        Note
        ----
        Allows quick subtraction of cubefiles (e.g. cfC=cfA+cfB)
        
        Parameters
        ----------
        other: cubefile
            cubefile to subtract
        
        Returns
        -------
        cubefile
            the subtraction 
        """
        return subtract(self, other)    
    
    def __mul__(self, factor):
        """
        Note
        ----
        Allows quick multiplication of a cubefiles and a float(e.g. cfC=cfA*factor)
        
        Parameters
        ----------
        factor: float
            the factor to multiply by
        
        Returns
        -------
        cubefile
            the product
        """
        return times(self, factor)    
     
    def __iadd__(self, other):
        """
        Note
        ----
        Allows addition of cubefiles (e.g. cfA+=cfB). NB: it changes cfA, returns nothing!
        
        Parameters
        ----------
        other: cubefile
            cubefile to add
        """
        self.add(other)
        return self    
    
    def __isub__(self, other):
        """
        Note
        ----
        Allows subtraction of cubefiles (e.g. cfA-=cfB). NB: it changes cfA, returns nothing!
        
        Parameters
        ----------
        other: cubefile
            cubefile to add
        """
        self.subtract(other)
        return self    

    def __imul__(self, factor):
        """
        Note
        ----
        Allows subtraction of cubefiles (e.g. cfA-=cfB). NB: it changes cfA, returns nothing!
        
        Parameters
        ----------
        factor: float
            multiplying factor
        """
        self.times(factor)
        return self  
        
    def check_same_grid(self, cf):
        """
        Parameters
        ----------
        cf: cubefile
            cubefile object to compare grid of
        """
        check_same_grid(self, cf)            
            
    def add(self, cf):
        """
        Parameters
        ----------
        cf: cubefile
            cubefile to add
        """
        self.check_same_grid(cf)
        self.values = np.add(self.values, cf.values)
        
    def subtract(self, cf):
        """
        Parameters
        ----------
        cf: cubefile
            cubefile to subtract
        """
        self.check_same_grid(cf)
        self.values = np.subtract(self.values, cf.values)
        
    def times(self, factor):
        """
        Parameters
        ----------
        factor: float
            factor to multiply by
        """
        self.values = np.multiply(self.values, factor)
        
    def trim(self, last=True):
        """
        Note
        ----
        Allows quick trimming of the last (default) or first value of a cube.
        
        Parameters
        ----------
        last: bool
            whether the last(True) or first (False) points of every axis should be trimmed.
        """
        if last == True:
            self.values = self.values[:-1,:-1,:-1]
        else:
            self.values = self.values[1:,1:,1:]
        self.Np_vect -= 1
        
    def RSR(self,cf):
        """       
        Parameters
        ----------
        cf: cubefile
            the cubefile to compare to
        
        Returns
        -------
        float
            the real space R value (RSR)
        """
        self.check_same_grid(cf)
        return RSR(self,cf)
    
#    def num(self,cf):
#        return num(self,cf)
#    def den(self,cf):
#        return den(self,cf)
        
    def get_coords(self, arr_pos): #nb. single point should be still given as array [[x],[y],[z]]
        """
        Note
        ----
        Allows to get all the coordinates as a grid. If single point, still provide as array [[x],[y],[z]]
        
        Parameters
        ----------
        arr_pos: array()
            the ensemble of points to get the coordinates of
        
        Returns
        -------
        array
            all the coordinates
        """
        all_coords = np.zeros(arr_pos.shape)
        for n,i in enumerate(arr_pos):
            coords = np.copy(self.origin)
            for j in range(3):
                coords += i[j]*self.Vect_M[j]
            all_coords[n] = coords
        return all_coords
    
    def get_pointlist(self):
        """
        Note
        ----print("Not the same grid! Exiting!")
            sys.exit()
        Allows to get the array of all array indexes.
        e.g.for a value array: np.arange(8).reshape(2,2,2)
        it would return:
            [[0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1]]
        
        Returns
        -------
        array
            the array of array indeces
        """
        return np.transpose(np.where(self.values != np.nan)) #trick to get it always verified
    
    def get_coordlist(self):
        """
        Returns
        -------
        the coordinates of the whole grid
        """
        return self.get_coords(self.get_pointlist())
    
    def min_dist(self, point):
        """
        Parameters
        ----------
        point: array(3,)
            point to measure distance from
        
        Returns
        -------
        float
            the shortest distance from the point to any atom of the geometry
        """
        d = np.zeros(len(self.geom.atoms))
        for n,i in enumerate(self.geom.coords):
            d[n] = np.linalg.norm(i - point)
        return d.min()
        
def read_cubefile(fname, empty_values=False):
    """
    Parameters
    ----------
    fname: str
        filename or path 
    empty_values: bool
        False => read all file, True => read only grid, zeros for values
    Returns
    -------
    cubefile 
        cubefile object
    """
    import mendeleev as md
    pnums = []
    V_m = []
    cmmnt = ""
    natm = 1        
    positions = []
    atms = []
    with open(fname) as f:
        for j,line in enumerate(f):
            if j < 2:
                cmmnt += line
            if j == 2:
                splt =line.split()
                natm = int(splt[0])
                o = np.asarray(list(map(float, splt[1:4])))
            if j > 2 and j < 6:
                splt =line.split()
                pnums.append(int(splt[0]))
                V_m.append(np.asarray(list(map(float, splt[1:4]))))
            if j >= 6 and j < natm+6:
                splt=line.split()
                atms.append(md.element(int(splt[0])).symbol)
                positions.append(np.asarray(list(map(float, splt[2:5]))))
            if j == natm+6:
                break
    Np_vect = np.asarray(list(map(int, pnums)))
    Vect_M = np.asarray(V_m)
    coords = np.asarray(positions)
    atoms = np.asarray(atms)
    g = geom(atoms, coords, coord_unit="au")
    if not empty_values:
        with open(fname) as f:
            cube = f.read()
        cubelist = []
        pattern2 = '(-?\d+\.?\d+[Ee]\ *[-+]?\ *\d+\s*){' + re.escape(str(pnums[2])) + '}'
        cubelist.append([])
        xcnt = 0
        ycnt = 0
        for xy,m in enumerate(re.finditer(pattern2,cube)):
            if xy % Np_vect[1] == 0 and xy > 0:
                cubelist[xcnt] = np.asarray(cubelist[xcnt])
                xcnt += 1 if xy > 0 else 0
                ycnt = 0
                cubelist.append([]) # append empty list for x
            s = m.group().split()
            s = np.asarray(list(map(float,s)))
            cubelist[xcnt].append(s)
            ycnt += 1
        cubelist[-1] = np.asarray(cubelist[-1])
        cubearray = np.asarray(cubelist)
    else:
        cubearray = np.zeros(list(Np_vect))
    return cubefile(g, Np_vect, Vect_M, cubearray, origin=o, comment=cmmnt)

def box_for_molecule(g="", margin=2.0):
    """
    Parameters
    ----------
    g: str or geom
        geometry (object of filename)
    margin: float or in
        the margin from the extreme nuclei
        
    Returns
    -------
    array(3,2)
        [[i.min,i.max] for i in x,y,z]
    """
    if g == "":
        g = input("g must be either geom object or .xyz filename. Type the filename, please\n")
    if type(g) == str:
        g = geom.from_xyz(g)
    c = g.copy()
    c.change_coord_unit("au")
    return np.array([c.coords.min(axis=0)-margin, c.coords.max(axis=0) + margin]).T

def grid_from_box(box, dist=0.2, fix="box"):
    """
    Parameters
    ----------
    box: array(3,2)
        [[i.min, i.max] for i in x,y,z] 
    dist: int or float
        distance between points
    fix: str
        "box" and some others keep box unchanged and slightly reduce the distance
        "dist" and some others keep the distance unchanged and slightly increase the box size
    
    Returns
    -------
    tuple(Np_vect, Vect_M, o)
        Vector of number of points, Vector size matrix, origin
    """
    box_size = box.ptp(axis=1)
    Np_vect = np.divide(box_size,dist)+1
    from math import ceil
    Np_vect = np.array([ceil(i) for i in Np_vect])
    o = box[:,0]
    if fix in ["box","extremes","margins"]:  # box is kept but the voxel is reduced slightly
        Vect_M = np.diag(np.divide(box_size, Np_vect))
        print("Your voxel now has size {},{},{}".format(*np.divide(box_size, Np_vect)))
    if fix in ["dist","distance", "voxel", "vect", "vector", "maxdist", "max_dist"]:
        Vect_M = np.diag(3*[dist])
        box[:,1] = box[:,0] + dist*Np_vect
        print("Your box has been changed to: ({},{}),({},{}),({},{})".format(*box.reshape(-1)))
    return (Np_vect, Vect_M, o)

def box_with_emptyval(g="", margin=2.0, dist=0.2, fix="box", comment="Empty function values\nOnly grid\n"):
    """
    Parameters
    ----------
    g: str or geom
        geometry (object of filename)
    margin: float or in
        the margin from the extreme nuclei
    dist: int or float
        distance between points
    fix: str
        "box" and some others keep box unchanged and slightly reduce the distance
        "dist" and some others keep the distance unchanged and slightly increase the box size    
    comment: str
        comment, default specifies it is empty
        
    Returns
    -------
    cubefile
        emptyval cubefile        
    """
    g = list(g) if type(g) == tuple else g  # no tuples because immutable!!
    if type(g) == str:
        if g == "":
            g = input("g must be either geom object or .xyz filename. Type the filename, please\n")
        g = [i.strip() for i in g.split(",")] if len(g.split(",")) == 2 else [g]  # Now g is a either [g1, g2] or [g]
    if type(g) == list:  # let's "convert" any filename to geom
        for n, i in enumerate(g):
            try:
                g[n] = geom.from_xyz(i, identifier=n)
            except:
                pass
        g = g[0] + g[1] if len(g) == 2 else g[0]
#    elif not isinstance(g, geom): 
    elif type(g) != geom: 
        raise TypeError("the geometry is neither a geometry object nor a file!")
    box = box_for_molecule(g, margin=margin)
    Np_vect, Vect_M, o =  grid_from_box(box, dist, fix=fix)
    return emptyval_cube_from_specs(g=g, Np_vect=Np_vect, Vect_M=Vect_M, o=o, comment=comment)
    
    
def emptyval_cube_from_specs(g="", Np_vect=np.array([]), Vect_M=np.array([]),  o=np.array([]), comment="Empty function values\nOnly grid\n"):
    """
    Note
    ----
    Allows to get a cubefile object with grid specs but empty values
    
    Parameters
    ----------
    g: geom obj,str, or list thereof
        geometry for the cubefile object or .xyz file to obtain it.
        if list, it combines the two geometries
    Np_vect: array(3,)
        number of points along x,y,z
    Vect_M: array(3,3)
        grid vectors. Generally diagonal
    o: array(3,)
        origin of the grid. default is 0,0,0
    comment: str
        comment, default specifies it is empty
        
    Returns
    -------
    cubefile
        emptyval cubefile
    """
    g = list(g) if type(g)==tuple else g  # no tuples because immutable!!
    if Np_vect == np.array([]):
        raise gridError("You must specify Np_vect")
    if Vect_M == np.array([]):
        raise gridError("You must specify Vect_M")
    if o == np.array([]):
        o = np.zeros(3)
    if type(g) == str:
        if g == "":
            g = input("g must be either geom object or .xyz filename. Type the filename, please\n")
        g = [i.strip() for i in g.split(",")] if len(g.split(",")) == 2 else [g]  # Now g is a either [g1, g2] or [g]
    if type(g) == list:  # let's "convert" any filename to geom
        for n,i in enumerate(g):
            try:
                g[n] = geom.from_xyz(i, identifier=n)
            except:
                pass
        g = g[0] + g[1] if len(g) == 2 else g[0]
#    elif not isinstance(g, geom): 
    elif type(g) != geom: 
        raise TypeError("the geometry is neither a geometry object nor a file!")
    values = np.zeros(list(Np_vect))
    g.change_coord_unit("au")
    return cubefile(g, Np_vect, Vect_M, values, o, comment=comment)
    
def cube_from_dm_and_specs(dmfile="", basfile="", g="", Np_vect=np.array([]), Vect_M=np.array([]),  o=np.array([]), comment=""):
    """
    Parameters
    ----------
    dmfile: str or list
        density matrix file, list or "file1,file2"
    basfile: str or list
        basis .nwchem file, list or "file1.nwchem,file2.nwchem"
    g: str or geom obj, or [str,str] or [obj,obj]
        geometry or geometry file which overrides self.geom. Compulsory only for dual basis. 
        For dual basis: "geomA.xyz, geomB.xyz",  ["geomA.xyz"," geomB.xyz"], or [geomA, geomB]
    Np_vect: array(3,)
        number of points along x,y,z
    Vect_M: array(3,3)
        grid vectors. Generally diagonal
    o: array(3,)
        origin of the grid. default is 0,0,0
        
    Returns
    -------
    cubefile object
        cube from the density matrix
    """
    cube = emptyval_cube_from_specs(g=g, Np_vect=Np_vect, Vect_M=Vect_M, o=o)
    cube.get_cubevals_from_dm(g=g, dmfile=dmfile, basfile=basfile)
    if comment:
        cube.comment = comment
    return cube

def cube_from_dm_and_box(dmfile="", basfile="", g="", margin=2.0, dist=0.2, fix="box", comment=""):
    """
    Parameters
    ----------
    dmfile: str or list
        density matrix file, list or "file1,file2"
    basfile: str or list
        basis .nwchem file, list or "file1.nwchem,file2.nwchem"
    g: str or geom obj, or [str,str] or [obj,obj]
        geometry or geometry file which overrides self.geom. Compulsory only for dual basis. 
        For dual basis: "geomA.xyz, geomB.xyz",  ["geomA.xyz"," geomB.xyz"], or [geomA, geomB]
    margin: float or in
        the margin from the extreme nuclei
    dist: int or float
        distance between points
    fix: str
        "box" and some others keep box unchanged and slightly reduce the distance
        "dist" and some others keep the distance unchanged and slightly increase the box size   
        
    Returns
    -------
    cubefile object
        cube from the density matrix
    """
    cube = box_with_emptyval(g=g, margin=margin, dist=dist, fix =fix)
    cube.get_cubevals_from_dm(g=g, dmfile=dmfile, basfile=basfile)
    if comment:
        cube.comment = comment
    return cube

def have_same_grid(cf1, cf2, orig_thresh=1e-5, printout=False):
    """
    Note
    ----
    Useful for scripts where different operations are to be done        
    
    Parameters
    ----------
    cf1: cubefile
        cubefile object 1
    cf2: cubefile
        cubefile object 2
    orig_thresh
        threshold to define the origins "equal", default is 1e-5
    printout: bool
        whether to print issue messages on screen
    Returns
    -------
    bool
        whether they have the same grid or not
    """
    if not (cf1.Np_vect == cf2.Np_vect).all():
        if printout:
            print("Not the same number of points per direction!!")
        return False
    elif not (cf1.Vect_M == cf2.Vect_M).all():
        if printout:
            print("Not the same grid vectors!!")
        return False
    elif not ((cf1.origin - cf2.origin) < orig_thresh).all():
        if printout:
            print("Not the same origin!!")
        return False
    elif cf1.values.shape != cf2.values.shape:
        if printout:
            print("Not the same number of total points!")
        return False
    else:
        return True
    
def check_same_grid(cf1, cf2, orig_thresh=0.00001):
    """
    Note
    ----
    Returns nothing! Just raises specific errors in case
    
    Parameters
    ----------
    cf1: cubefile
        cubefile object 1
    cf2: cubefile
        cubefile object 2
    orig_thresh
        threshold to define the origins "equal", default is 1e-5
    printout: bool
        whether to print issue messages on screen
    Returns
    -------
    bool
        whether they have the same grid or not
    """
    if not (cf1.Np_vect == cf2.Np_vect).all():
        raise gridError("number of points does not match")
    elif not (cf1.Vect_M == cf2.Vect_M).all():
        raise gridError("vectors do not match")
    elif not ((cf1.origin - cf2.origin) < orig_thresh).all():
        raise gridError("origins do not match")
    elif cf1.values.shape != cf2.values.shape:
        raise gridError("Number of values does not match")
        
def subtract(cf1, cf2, comment="subtracted with cubetools\n here it goes\n", geom="", Continue=False):
    """
    Parameters
    ----------
    cf1: cubefile1
        cubefile to subtract from
    cf2: cubefile2
        cubefile to subtract
    comment: str
        comment for the resulting cubefile object. Must be less than 2 lines
    geom: geom
        geometry object. Default uses cf1.geom
    Continue: bool
        whether to raise gridError or create mock-object in case of different grid.
        default is False, True is useful for "for" loops
    
    Returns
    -------
    cubefile
        cubefile object with difference
    """
    if comment != "subtracted with cubetools\n here it goes\n":
        if len(comment.split("\n")) == 3:
            pass
        elif len(comment.split("\n")) < 3:
            comment += "\n"*(3-len(comment.split("\n")))
        elif len(comment.split("\n")) > 3:
            raise otherError("Please retry with a two- or fewer-lines comment")
    if geom == "":
        geom = cf1.geom
    if have_same_grid(cf1, cf2):
        Np_vect = cf1.Np_vect
        Vect_M = cf1.Vect_M
        origin = cf1.origin
    else:
        if Continue:
            print("Not the same grid! creating mock-object and continuing")
            Np_vect = cf1.Np_vect
            Vect_M = cf1.Vect_M
            origin = cf1.origin
            return cubefile(geom, Np_vect, Vect_M, np.zeros(cf1.values.shape), origin, comment=comment)
        else:
            raise gridError()
    return cubefile(geom, Np_vect, Vect_M, np.subtract(cf1.values, cf2.values), origin, comment=comment)

def add(cf1, cf2, comment="added with cubetools\n here it goes\n", geom="", Continue=False):
    """
    Parameters
    ----------
    cf1: cubefile1
        cubefile to add to
    cf2: cubefile2
        cubefile to add
    comment: str
        comment for the resulting cubefile object. Must be less than 2 lines
    geom: geom
        geometry object. Default uses cf1.geom
    Continue: bool
        whether to raise gridError or create mock-object in case of different grid.
        default is False, True is useful for "for" loops
    
    Returns
    -------
    cubefile
        cubefile object with addition
    """
    if comment != "added with cubetools\n here it goes\n":
        if len(comment.split("\n")) == 3:
            pass
        elif len(comment.split("\n")) < 3:
            comment += "\n"*(3-len(comment.split("\n")))
        elif len(comment.split("\n")) > 3:
            raise otherError("Please retry with a two- or fewer-lines comment")
    if geom == "":
        geom = cf1.geom
    if have_same_grid(cf1, cf2):
        Np_vect = cf1.Np_vect
        Vect_M = cf1.Vect_M
        origin = cf1.origin
    else:
        if Continue:
            print("Not the same grid! creating mock-object and continuing")
            Np_vect = cf1.Np_vect
            Vect_M = cf1.Vect_M
            origin = cf1.origin
            return cubefile(geom, Np_vect, Vect_M, np.zeros(cf1.values.shape), origin, comment=comment)
        else:
            raise gridError()
    return cubefile(geom, Np_vect, Vect_M, np.add(cf1.values, cf2.values), origin, comment=comment)

def times(cf, factor):
    """
    Parameters
    ----------
    cf: cubefile
        cubefile to multiply
    factor: float
        factor to multiply by
    
    Returns
    -------
    cubefile
        cubefile object where values have been multiplied
    """
    copy = cf.copy()
    copy.times(factor)
    return copy
        
def RSR(cf1, cf2, Continue=False):
    """       
    Parameters
    ----------
    cf1: cubefile
        the cubefile to compare
    cf2: cubefile
        the cubefile to compare to
    
    Returns
    -------
    float
        the real space R value (RSR)
    """
    if not check_same_grid(cf1, cf2):
        if Continue:
            print("Not the same grid, but continuing!Returning 2!")
            return 2.0
        else:
            raise gridError()
    else:            
        num=(abs(np.subtract(cf1.values, cf2.values))).sum()
        den=(abs(np.add(cf1.values, cf2.values))).sum()
    return np.divide(num, den)


#def num(cf1,cf2,Continue=False):
#    if not check_same_grid(cf1,cf2):
#        if Continue:
#            print("Not the same grid, but continuing!Returning 2!")
#            return 2.0
#        else:
#            print("Not the same grid! Exiting!")
#            sys.exit()
#    else:            
#        num=(abs(np.subtract(cf1.values,cf2.values))).sum()
#    return num
#
#def den(cf1,cf2,Continue=False):
#    if not check_same_grid(cf1,cf2):
#        if Continue:
#            print("Not the same grid, but continuing!Returning 2!")
#            return 2.0
#        else:
#            print("Not the same grid! Exiting!")
#            sys.exit()
#    else:            
#        den=(abs(np.add(cf1.values,cf2.values))).sum()
#    return den

    


    