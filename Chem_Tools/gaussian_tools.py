
import sys
sys.path.append("/Users/coire/McCoy/") #Path to McUtils, References..etc
import numpy as np
import copy
import McUtils.GaussianInterface as GI
import os
import matplotlib.pyplot as plt
from Chem_Tools import misc, density_cube as dc


class GaussianResults:
    '''
    The purpose of this class is to combine the results from fchk and log files into a single object
    What gets pulled from where depends on the type of job run
    I am assuming this class will be updated a lot moving forward
    '''
    def __init__(self, files, **kwargs):
        self.params = kwargs
        if len(files) == 0:
            raise Exception('Nothing to interpret.')
        self.files = files
        if not isinstance(files, list): self.files = [self.files]
        self._MP2Energy = None
        self._input_zmat = None
        self.force_constants = None

    @property
    def MP2Energy(self):
        if self._MP2Energy is None:
            self._MP2Energy = self.get_MP2energy()
        return self._MP2Energy

    @property
    def input_zmat(self):
        if self._input_zmat is None:
            self._input_zmat = self.get_input_zmat()
        return self._input_zmat

    def get_MP2energy(self):
        energy = []
        for file in self.files:
            with GI.GaussianFChkReader(file + '.fchk') as reader:
                parse = reader.parse("MP2 Energy")
            ens = parse["MP2 Energy"]
            energy.append(ens)
        e = np.array(energy)
        return e

    def get_input_zmat(self):
        zmats = []
        for f in self.files:
            with GI.GaussianLogReader(f + '.log') as reader:
                parse = reader.parse("InputZMatrix")
            if parse["InputZMatrix"] is None:
                print(f'could not pull the z-matrix from {f}')
                print(parse["InputZMatrix"])
            new_zmat = misc.zmat(g_str=parse["InputZMatrix"])
            zmats.append(new_zmat)
        return zmats

    def pull_forces(self):
        forces = []
        for f in self.files:
            with GI.GaussianLogReader(f + '.log') as reader:
                parse = reader.parse("ForceConstants")
            forces.append(parse["ForceConstants"])
        self.force_constants=forces



class GaussianJob:
    def __init__(self, params, filename):
        self.params = params
        self.filename=filename

    def save_job(self):
        with open(self.filename, 'w', newline='\n') as f:
            for key in self.params:
                f.write(self.params[key])






class GLogInterpreter:
    def __init__(self, log_files):
        self.log_files = log_files
        if not isinstance(log_files, list): self.log_files = [self.log_files]
        self._Energy = None
        self._intCoords = None
        self._dipole_moms = None
        self._input_zmat = None
        self.force_constants = None
        self.NBO_Summary = None

    @property
    def Energy(self):
        if self._Energy is None:
            self._Energy = self.get_Ens()
        return self._Energy

    @property
    def intCoords(self):
        if self._intCoords is None:
            self._intCoords = self.get_intCoords()
        return self._intCoords

    @property
    def dipole_moms(self):
        if self._dipole_moms is None:
            self._dipole_moms = self.get_dipole_moms()
        return self._dipole_moms

    @property
    def input_zmat(self):
        if self._input_zmat is None:
            self._input_zmat = self.get_input_zmat()
        return self._input_zmat

    def get_Ens(self):
        energy_array = []
        for f in self.log_files:
            with GI.GaussianLogReader(f) as reader:
                parse = reader.parse("OptimizedScanEnergies")
            temp_energy_array, temp_coords = parse["OptimizedScanEnergies"]
            if np.any(energy_array):
                energy_array = np.append(energy_array, temp_energy_array)
            else:
                energy_array = copy.deepcopy(temp_energy_array)
        return energy_array

    def get_intCoords(self):
        coords = {}
        for f in self.log_files:
            with GI.GaussianLogReader(f) as reader:
                parse = reader.parse("OptimizedScanEnergies")
            temp_energy_array, temp_coords = parse["OptimizedScanEnergies"]
            if coords:
                for key in coords.keys():
                    coords[key] = np.append(coords[key],temp_coords[key])
            else:
                coords = copy.deepcopy(temp_coords)
        return coords

    def get_dipole_moms(self):
        dipoles = {}
        for f in self.log_files:
            with GI.GaussianLogReader(f) as reader:
                parse = reader.parse("DipoleMoments")
            #print(parse["DipoleMoments"])
            temp_dipoles = parse["DipoleMoments"]

            if dipoles:
                for key in dipoles.keys():
                    dipoles[key] = np.append(dipoles[key],temp_dipoles[key])
            else:
                dipoles = copy.deepcopy(temp_dipoles)
        return dipoles

    def get_input_zmat(self):
        zmats = []
        for f in self.log_files:
            with GI.GaussianLogReader(f) as reader:
                parse = reader.parse("InputZMatrix")
            zmats.append(misc.zmat(g_str= parse["InputZMatrix"]))
        return zmats

    def pull_forces(self):
        for f in self.log_files:
            with GI.GaussianLogReader(f) as reader:
                parse = reader.parse("ForceConstants")
            forces= parse["ForceConstants"]
        self.force_constants=forces

    def pull_NBO(self):
        d = {}
        ex = {'occupancy': [],
                 'file': []}
        for f in self.log_files:
            with GI.GaussianLogReader(f) as reader:
                parse = reader.parse("NBOSummary")
            NBO= parse["NBOSummary"]
            for o in range(len(NBO[0])):
                if NBO[0][o] in d:
                    d[NBO[0][o]]['occupancy'].append(NBO[1][o])
                    d[NBO[0][o]]['file'].append(f)
                else:
                    new_d = copy.deepcopy(ex)
                    new_d['occupancy'].append(NBO[1][o])
                    new_d['file'].append(f)
                    d[NBO[0][o]] = new_d
        self.NBO_Summary = d





class GLogPlotter:
    def __init__(self, files):
        self.files = files
        self._data = None

    @property
    def data(self):
        if self._data is None:
            self._data = self.get_data()
        return self._data

    def get_data(self):
        data = GLogInterpreter(self.files)
        return data

    def plot_energies(self):
        x = self.data.Energy
        fig = plt.figure()
        plt.plot(x, '.')

    def plot_energy_v_coord(self, xCoord):
        y = self.data.Energy
        x = self.data.intCoords[xCoord]
        plt.plot(x, y, ".")

    def plot_energy_v_twoCoords(self, xCoord, yCoord):
        x = self.data.intCoords[xCoord]
        y = self.data.intCoords[yCoord]
        z = self.data.Energy
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter3D(x, y, z, 'gray')
        ax.set_xlabel(xCoord)
        ax.set_ylabel(yCoord)
        ax.set_zlabel('Energy')

    def plot_coord_v_coord(self, xCoord, yCoord):
        x = self.data.intCoords[xCoord]
        y = self.data.intCoords[yCoord]
        fig = plt.figure()
        ax = plt.axes()
        ax.plot(x,y,'.')
        ax.set_xlabel(xCoord)
        ax.set_ylabel(yCoord)

class FchkInterpreter:
    # Copied directly from https://github.com/rmhuch/QOOH/blob/master/FChkInterpreter.py
    #then modified
    def __init__(self, fchks, cubes=False,  **kwargs):
        self.params = kwargs
        if len(fchks) == 0:
            raise Exception('Nothing to interpret.')
        self.fchks = fchks
        if not isinstance(fchks, list): self.fchks = [self.fchks]
        self._hessian = None
        self._cartesians = None  # dictionary of cartesian coordinates keyed by (x, y) distances
        self._gradient = None
        self._MP2Energy = None
        self._atomicmasses = None
        self._atomicNumbers = None
        self._mullikens = None
        self.cubes= cubes
        if self.cubes:
            self.initializeCubes()

    @property
    def cartesians(self):
        if self._cartesians is None:
            self._cartesians = self.get_coords()
        return self._cartesians

    @property
    def hessian(self):
        if self._hessian is None:
            self._hessian = self.get_hess()
        return self._hessian

    @property
    def gradient(self):
        if self._gradient is None:
            self._gradient = self.get_grad()
        return self._gradient

    @property
    def MP2Energy(self):
        if self._MP2Energy is None:
            self._MP2Energy = self.get_MP2energy()
        return self._MP2Energy

    @property
    def atomicmasses(self):
        if self._atomicmasses is None:
            self._atomicmasses = self.get_mass()
        return self._atomicmasses

    @property
    def atomicNumbers(self):
        if self._atomicNumbers is None:
            self._atomicNumbers = self.get_atomicNumbers()
        return self._atomicNumbers

    @property
    def mullikens(self):
        if self._mullikens is None:
            self._mullikens = self.get_mullikens()
        return self._mullikens

    def get_coords(self):
        """Uses McUtils parser to pull cartesian coordinates
            :returns coords: nx3 coordinate matrix"""
        crds = []
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("Coordinates")
            coords = parse["Coordinates"]
            crds.append(coords)
        c = np.array(crds)
        if c.shape[0] == 1:
            c = np.squeeze(c)
        return c

    def get_hess(self):
        """Pulls the Hessian (Force Constants) from a Gaussian Frequency output file
            :arg chk_file: a Gaussian Frequency formatted checkpoint file
            :returns hess: full Hessian of system as an np.array"""
        hess = []
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("ForceConstants")
            forcies = parse["ForceConstants"]
            hess.append(forcies)
        h = np.array(hess)
        # ArrayPlot(forcies.array, colorbar=True).show()
        return h

    def get_grad(self):
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("Gradient")
            grad = parse["Gradient"]
        return grad

    def get_MP2energy(self):
        energy = []
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("MP2 Energy")
            ens = parse["MP2 Energy"]
            energy.append(ens)
        e = np.array(energy)
        return e

    def get_mass(self):
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("AtomicMasses")
            mass_array = parse["AtomicMasses"]
        return mass_array

    def get_atomicNumbers(self):
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("AtomicNumbers")
            nums = parse["AtomicNumbers"]
        return nums

    def get_mullikens(self):
        mulls = []
        for fchk in self.fchks:
            with GI.GaussianFChkReader(fchk) as reader:
                parse = reader.parse("Mulliken Charges")
            mulls.append(parse["Mulliken Charges"])
        return np.array(mulls)

    def initializeCubes(self):
        ed_cds = []
        ed_vals = []
        for f in self.fchks:
            print(f)
            cube = dc.cube(os.path.splitext(f)[0] + '.cube')
            ed_cds.append(cube.cds)
            ed_vals.append(cube.density)
        self.ed_cds = np.array(ed_cds)
        self.ed_vals = np.array(ed_vals)



#class NormalModes:
    #def __init__(self):

#class ZMatrix:

class CoirePlots:
    def __init__(self, *data, type=None):
        self.data = data
        self.type = type
        self._axSettings = None
        self._pltSettings = None

    @property
    def axSettings(self):
        if self._axSettings is None:
            self._axSettings = self.emptyAxSettings()
        return self._axSettings

    @property
    def pltSettings(self):
        if self._pltSettings is None:
            self._pltSettings = self.emptyPltSettings()
        return self._pltSettings

    def emptyAxSettings(self):
        settings = {'xlabel': None,
                    'ylabel': None,
                    'title': None}
        return settings

    def plotType(self, type):
        self.type = type

    #def setSettings(self, ax):

    def plot(self, **kwargs):
        if self.type == 'scatter2d':
            self.fig = plt.figure()
            self.ax = plt.axes()
            self.ax.scatter(self.data[0], self.data[1], kwargs)
            if self.axSettings is not None:
                self.setAxSettings()
            if self.pltSettings is not None:
                self.setPltSettings()

    def setAxSettings(self):
        self.ax.set_xlabel(self.axSettings['xlabel'])
        self.ax.set_ylabel(self.axSettings['xlabel'])
        self.ax.title(self.axSettings['title'])


