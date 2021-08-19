import numpy as np
from scipy import linalg
import utilities as uts
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import energy_arrays
import math

"""
Future DVR implementations (nD DVRs) will be placed here
"""
class DVR_1D:
    def __init__(self,
                 grid,
                 kinetic_matrix,
                 potential_matrix,
                 filename,
                 save_method='pickle'):
        """
        M_1D based on Colbert-Miller DVR described here: http://www.cchem.berkeley.edu/millergrp/pdf/243.pdf
        To use this DVR, supply the class with a grid, a kinetic energy matrix, and a potential energy
        matrix. The DVR does the incredibly complex step of saying H = T + V and diagonalizing the resulting H matrix.
        The actually thinking happens when constructing the matrices...
        The results are put in a dictionary and saved by the method of your choosing (pickling or npz, right now just pickle)

        :param grid: the CM 1d DVR grid. and equally-spaced array of points along your coordinate of interest
        :type grid: np.ndarray
        :param kinetic_matrix: the matrix representation of your kinetic energy operator in the CM-DVR basis
        :type kinetic_matrix: np.ndarray
        :param potential_matrix: the matrix representation of your potential energy. This is the electronic potential energy of your system at each grid point.
        :type potential_matrix: np.ndarray
        :param filename: what to name the results file
        :type filename: str
        :param save_method: how to save the results file (e.g. 'pickle', 'npz', etc.)
        :type save_method: str
        """
        self.grid = grid
        self.kinetic_matrix = kinetic_matrix
        self.potential_matrix = potential_matrix
        self.filename = filename
        self.save_method = save_method

    def run(self):
        """
        H = T + V
        returns the eigenvalues and eigenvectors of the diagonalized matrix as energies and wfns in a dictionary
        """
        hamiltonian_matrix = self.kinetic_matrix + self.potential_matrix
        eigvals, eigvecs = linalg.eigh(hamiltonian_matrix)
        results = {"grid": self.grid,
                   "potential":np.diag(self.potential_matrix),
                   "kinetic":self.kinetic_matrix,
                   "energies": eigvals,
                   "wfns": eigvecs}
        uts.save(object=results, method= self.save_method, filename=self.filename)

class Grid:
    def __init__(self,
                 domain,
                 num_points,
                 inclusive=True,):
        """
        makes an equally-spaced grid on the interval [domain[0],domain[1]] if inclusive=True or [domain[0],domain[1]) if inclusive=False
        :param domain: [start, stop]
        :type domain: tuple
        :param num_points: the number of grid points that you want
        :type num_points: int
        :param inclusive: do you want to include the last point in the grid?
        :type inclusive: bool
        """
        self.domain = domain
        self.num_points = num_points
        self.inclusive = inclusive
        self.initialize()

    def initialize(self):
        """
        makes the grid using np.linspace
        np.linspace has the native argument "endpoint" which include the endpoint in the resulting grid only if true
        """
        self.grid = np.linspace(self.domain[0], self.domain[1], self.num_points, endpoint=self.inclusive)

class AnalyzeDVR:
    def __init__(self,
                 results_file,
                 save_method='pickle',
                 in_AU = True,
                 energy_units='wavenumbers',
                 grid_unit='angstroms'):
        """
        reads DVR results from a file for analysis. All unit conversions are done when the file is read and all subsequent analysis methods assume the data is in the form that you want.
        :param results_file: the file that has the DVR results
        :type results_file: str
        :param save_method: the method used to save the DVR results (e.g. 'pickle, 'npz', etc.)
        :type save_method: str
        :param in_AU: Are the DVR results loaded in atomic units, default is True bc my DVR code only deals with atomic units
        :type in_AU: bool
        :param energy_units: the energy units you want to deal with. Default is 'wavenumbers'
        :type energy_units: str
        :param grid_unit: the units you want your DVR grid in (e.g. 'angstroms', 'radians', etc).
        """
        self.results_file = results_file
        self.save_method = save_method
        self._results = None
        self.in_AU = in_AU
        self.energy_units=energy_units
        self.grid_unit=grid_unit
        if in_AU:
            self.convert_units()

    @property
    def results(self):
        if self._results is None:
            self._results = self.get_results()
        return self._results

    def get_results(self):
        results = uts.load(filename=self.results_file,method=self.save_method)
        for i in range(len(results['wfns'])):
            wfn = results['wfns'][:,i]
            phase = np.sign(wfn[np.argmax(np.abs(wfn))])
            results['wfns'][:,i] = results['wfns'][:,i] * phase
        return results

    def get_frequency(self, excitation):
        freq_au = self.results['energies'][excitation[1]] - self.results['energies'][excitation[0]]
        return freq_au

    def convert_units(self):
        if self.grid_unit != 'radians':
            self.results['grid'] = uts.Constants.convert(self.results['grid'], self.grid_unit, to_AU=False)
        self.results['energies'] = uts.Constants.convert(self.results['energies'], self.energy_units, to_AU=False)
        self.results['potential'] = uts.Constants.convert(self.results['potential'], self.energy_units, to_AU=False)
        self.results['kinetic'] = uts.Constants.convert(self.results['kinetic'], self.energy_units, to_AU=False)


    def plot_wfns(self, states, on_pot=False, x_range=[], y_range=[], scale=1, save_file=''):
        if not len(x_range) > 1:
            x_range = [self.results['grid'].min(), self.results['grid'].max()]
        x_argmin = np.argmax(self.results['grid'] >= x_range[0])
        x_argmax = np.argmin(self.results['grid'] <= x_range[1]) - 1
        if not isinstance(states, list): states = [states]

        evenly_spaced_interval = np.linspace(0,0.9, len(states))
        colors = [cm.viridis(x) for x in evenly_spaced_interval]

        fig = plt.figure()
        ax = plt.axes()
        i=0
        for state in states:
            ax.plot(self.results['grid'][x_argmin:x_argmax],
                    (self.results['wfns'][:,state][x_argmin:x_argmax] * scale) + self.results['energies'][state],
                    label=r'$\psi_{%s}$' %state,
                    color=colors[i])
            ax.axhline(y=self.results['energies'][state], color=colors[i], linestyle='--')
            i += 1

        if on_pot:
            ax.plot(self.results['grid'][x_argmin:x_argmax],
                    self.results['potential'][x_argmin:x_argmax],
                    label="potential", color='black',
                    linewidth=2,
                    linestyle='--')

        ax.set_xlim(x_range)
        if len(y_range) > 0:
            ax.set_ylim(y_range)
        #ax.set_xlabel(r'$\Delta r$')
        plt.xlabel("rxn coordinate (Angstroms)")
        plt.ylabel("Energy (cm-1)", labelpad=15)
        plt.legend()
        plt.tight_layout()
        if save_file:
            plt.savefig(save_file)
        plt.show()

    def plot_energies(self, energies_range=[], energy_step=1, save_file='', x_range=[], y_range=[]):
        if not len(x_range) > 1:
            x_range = [self.results['grid'].min(), self.results['grid'].max()]
        x_argmin = np.argmax(self.results['grid'] >= x_range[0])
        x_argmax = np.argmin(self.results['grid'] <= x_range[1]) - 1

        if not len(energies_range) > 0:
            energies_range = [0, len(self.results['energies'])]
        energies_to_plot = np.arange(energies_range[0], energies_range[1], energy_step)

        evenly_spaced_interval = np.linspace(0,0.9, len(energies_to_plot))
        colors = [cm.viridis(x) for x in evenly_spaced_interval]

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(self.results['grid'][x_argmin:x_argmax],
                self.results['potential'][x_argmin:x_argmax],
                label="potential", color='black',
                linewidth=2,
                linestyle='--')
        i=0
        for e in energies_to_plot:
            ax.axhline(y=self.results['energies'][e],
                       color=colors[i])
            i+=1
        ax.set_xlim(x_range)
        if len(y_range) > 0:
            ax.set_ylim(y_range)
        plt.xlabel("rxn coordinate (Angstroms)")
        plt.ylabel("Energy (cm-1)", labelpad=15)
        plt.tight_layout()
        if save_file:
            plt.savefig(save_file)
        plt.show()

    def plot_ind_wfns(self, wfns_range=[0,8], wfns_step=1, scale=1, x_range=[], y_range=[], save_file=''):
        if not len(x_range) > 1:
            x_range = [self.results['grid'].min(), self.results['grid'].max()]
        x_argmin = np.argmax(self.results['grid'] >= x_range[0])
        x_argmax = np.argmin(self.results['grid'] <= x_range[1]) - 1
        wfns_to_plot = np.arange(wfns_range[0], wfns_range[1]+1, wfns_step)
        rows = math.ceil(len(wfns_to_plot)/3)

        evenly_spaced_interval = np.linspace(0, 0.9, len(wfns_to_plot))
        colors = [cm.viridis(x) for x in evenly_spaced_interval]

        fig, axs = plt.subplots(rows, 3, sharey=True, sharex=True)
        i=0
        for r in range(rows):
            for c in range(3):
                if i < len(wfns_to_plot):
                    state = wfns_to_plot[i]
                    axs[r,c].plot(self.results['grid'][x_argmin:x_argmax],
                            self.results['potential'][x_argmin:x_argmax],
                            color='black',
                            linewidth=2,
                            linestyle='--')
                    axs[r,c].plot(self.results['grid'][x_argmin:x_argmax],
                                  (self.results['wfns'][:, state][x_argmin:x_argmax] * scale) + self.results['energies'][state],
                                  label=r'$\rm{\psi_{%s}}$' % state,
                                  color=colors[i])
                    axs[r,c].axhline(y=self.results['energies'][state], color=colors[i], linestyle='--')
                    axs[r,c].set_xlim(x_range)
                    axs[r,c].legend()
                    axs[r,c].set_xlim(x_range)
                    if len(y_range) > 0:
                        axs[r,c].set_ylim(y_range)
                    i+=1
        plt.setp(axs, ylim=axs[rows-1, 2].get_ylim())
        # add a big axis, hide frame
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel("rxn coordinate (Angstroms)")
        plt.ylabel("Energy (cm-1)", labelpad=15)
        #plt.subplots_adjust(wspace=0, hspace=0)
        plt.tight_layout()
        if save_file:
            plt.savefig(save_file)
        plt.show()

    def get_zpe(self):
        return self.results['energies'][0]






#if __name__ == "__main__":


    """OH_mass = ut.Constants.reduced_mass('O-H', to_AU=True)
    grid_range = ut.Constants.convert(1, 'angstroms', to_AU=True)
    freq = ut.Constants.convert(3300, 'wavenumbers', to_AU=True)

    grid = Grid((-grid_range,grid_range), 100).grid
    kin_mat = energy_arrays.CM_1D_kin_mat(grid, interval='-infty_to_infty', mass=OH_mass)
    kinetic_matrix = kin_mat.matrix"""

    """pot_opts = {'type': 'harmonic oscillator',
                'mass': OH_mass,
                'frequency': freq}

    pot_mat = energy_arrays.PredefinedPotMats(grid, pot_opts)
    potential_matrix = pot_mat.matrix
    """
    """gaussian_grid = np.flip(np.array([1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5])) -0.95785
    gaussian_grid = ut.Constants.convert(gaussian_grid, "angstroms", to_AU=True)
    gaussian_energies = np.flip(np.array([-0.30799, -0.32956, -0.35195, -0.37395, -0.39332, -0.40575,
                                          -0.40329, -0.37102, -0.28049, -0.07581, 0.35845])) + -76.0
    scan_results = np.array([gaussian_grid, gaussian_energies])


    pot_mat = energy_arrays.ScanPotMat(grid=grid, scan_results=scan_results, min_shift=True)
    potential_matrix = pot_mat.matrix"""
    """CH3_mass = uts.Constants.reduced_mass('C-H-H-H', to_AU=True)
    grid = Grid((0, 2 * np.pi), 101, inclusive=False).grid

    kin_mat = energy_arrays.CM_1D_kin_mat(grid, interval='0_to_2pi', mass=0.5)
    kinetic_matrix = kin_mat.matrix
    #plt.matshow(kinetic_matrix)

    pot_opts = {'type': 'periodic',
                'alpha':uts.Constants.convert(200, 'wavenumbers', to_AU=True)}
    pot_mat = energy_arrays.PredefinedPotMats(grid=grid, pot_opts=pot_opts)
    #plt.plot(grid, pot_mat.array)
    #potential_matrix = pot_mat.matrix
    potential_matrix = np.diag(np.zeros(len(grid)))

    file = '/Users/coire/McCoy/QOOH_repo/coire/DVR_dev/water_monomer_rpath/DVR_test'
    test = DVR_1D(grid=grid, kinetic_matrix=kinetic_matrix, potential_matrix=potential_matrix, filename=file)
    test.run()

    anl_test = AnalyzeDVR(results_file=file, in_AU=False, energy_units='wavenumbers', grid_unit='radians')
    anl_test.get_frequency((0,1))
    anl_test.rephase_wfns([0])
    anl_test.plot_wfns(states=[0, 1, 3],
                       pot=np.diag(potential_matrix),
                       on_pot=False)
                       #x_range=[-0.4, 0.5],
                       #y_range=[-2500,12000],
                       #scale=0.03,
                       #save_file=file + '_0_2pi.png')
    zpe = anl_test.get_zpe()
    test = DVR_1D(grid=grid, kinetic_matrix=kinetic_matrix, potential_matrix=potential_matrix, filename=file)
    test.run()
    print("hello")"""



