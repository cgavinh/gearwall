import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as la
import sys
sys.path.append("../../")
sys.path.append("..")



class cube:
    def __init__(self,
                 filename):
        self._cube = None
        self.filename = filename
        self.initialize()

    def initialize(self):
        atmDict = {1: 'H',
                   3: 'Li',
                   8: 'O',
                   11: 'Na',
                   55: 'Cs',
                   17: 'Cl',
                   35: 'Br',
                   53: 'I'}
        a = open(self.filename, "r")
        lin = a.readlines()
        beginningE = int(lin[2].split()[0]) + 6
        linNum = 1
        atmCt = 0
        self.atmStr = []
        ndim = []
        delta = []
        gridN = 0
        for line in lin:
            if linNum == 3:
                splt = line.split()
                info = [float(q) for q in splt]
                nAtoms = int(info[0])
                self.atmAr = np.zeros((nAtoms, 3))
                # atmStr = np.tile(["O","H","H"],nAtoms//3)
                origin = np.array(info[1:-1])
            elif 4 <= linNum <= 6:
                splt = line.split()
                xyzZ = [float(q) for q in splt]
                ndim.append(int(xyzZ[0]))
                if linNum == 4:
                    delta = xyzZ[1]
            elif 7 <= linNum <= beginningE:
                splt = line.split()
                xyzZ = [float(q) for q in splt][1:]
                self.atmStr.append(atmDict[int(splt[0])])
                self.atmAr[atmCt] = xyzZ[1:]
                atmCt += 1
            elif linNum > beginningE:
                pass
            linNum += 1
        a.close()
        k = pd.read_table(self.filename, delim_whitespace=True, header=None, skiprows=np.arange(beginningE))
        xx = k.to_numpy()
        xxp = xx[np.logical_not(np.isnan(xx))]
        self.density = np.reshape(xxp, ndim)
        self.cdsX = np.linspace(origin[0], origin[0] + delta * ndim[0], num=ndim[0])
        self.cdsY = np.linspace(origin[1], origin[1] + delta * ndim[1], num=ndim[1])
        self.cdsZ = np.linspace(origin[2], origin[2] + delta * ndim[2], num=ndim[2])
        self.cds = {'x' : self.cdsX, 'y': self.cdsY, 'z':self.cdsZ}
        self.arr = np.meshgrid(self.cdsX, self.cdsY, self.cdsZ, indexing='ij')
        dxO = self.cdsX[1] - self.cdsX[0]
        dyO = self.cdsY[1] - self.cdsY[0]
        dzO = self.cdsZ[1] - self.cdsZ[0]
        print("DXDYDZO", dxO, dyO, dzO, dxO * dyO * dzO)
        nbins = 40
        shiftL = 1.0
        shiftR = 1.0
        cylRad = 1.0
