import numpy as np
from scipy import linalg
import utilities as ut

#CM_1D based on Colbert-Miller DVR described here: http://www.cchem.berkeley.edu/millergrp/pdf/243.pdf
class DVR_1D:
    def __init__(self,
                 grid,
                 kinetic_matrix,
                 potential_matrix,
                 file_name,
                 save_method='pickle'):

        self.grid = grid
        self.kinetic_matrix = kinetic_matrix
        self.potential_matrix = potential_matrix
        self.file_name = file_name
        self.save_method = save_method

    def run(self):
        hamiltonian_matrix = self.kinetic_matrix + self.potential_matrix
        energies, wfn_coeffs = linalg.eigh(hamiltonian_matrix)
        results = {"grid": self.grid, "energies": energies, "wfns": wfn_coeffs}
        self.save(results)

    def save(self, results):
        if self.save_method == 'pickle':
            ut.pickle_save(results, self.file_name)

class Grid:
    def __init__(self,
                 domain,
                 num_points,
                 inclusive=True):
        self.domain = domain
        self.num_points = num_points
        self.inclusive = inclusive
        self.initialize()

    def initialize(self):
        """
        np.linspace has input endpoint, if True, it includes stop, else excludes it
        :return:
        """
        self.grid = np.linspace(self.domain[0], self.domain[1], self.num_points, endpoint=self.inclusive)

class AnalyzeDVR:
    def __init__(self,
                 results_file,
                 file_type='pickle',
                 in_AU = True,
                 energy_units='wavenumbers',
                 grid_unit='angstroms'):
        self.results_file = results_file
        self.file_type = file_type
        self._results = None
        self.in_AU = in_AU
        self.energy_units=energy_units
        self.grid_unit=grid_unit
        if not in_AU:
            self.convert_units()

    @property
    def results(self):
        if self._results is None:
            self._results = self.get_results()
        return self._results

    def get_results(self):
        results = ut.pickle_load(self.results_file)
        return results

    def get_frequency(self, excitation):
        freq_au = self.results['energies'][excitation[1]] - self.results['energies'][excitation[0]]
        return freq_au

    def convert_units(self):
        if self.grid_unit != 'radians':
            self.results['grid'] = ut.Constants.convert(self.results['grid'], self.grid_unit, to_AU=False)
            self.results['energies'] = ut.Constants.convert(self.results['energies'], self.energy_units, to_AU=False)
        self.results['wfns'] = ut.Constants.convert(self.results['wfns'], self.energy_units, to_AU=False)

    def rephase_wfns(self, states):
        """
        states is a list of the wfns that need to be re-phased
        :param states:
        :return:
        """
        phases = np.ones(len(self.results['grid']))
        for state in states:
            phases[state] = -1
        self.results['wfns'] = self.results['wfns']*phases

    def plot_wfns(self, states, pot = [], on_pot=False, x_range=[], y_range=[], scale=1, save_file=''):
        if not len(x_range) > 1:
            x_range = [self.results['grid'].min(), self.results['grid'].max()]
        x_min = np.argmax(self.results['grid'] >= x_range[0])
        x_max = np.argmin(self.results['grid'] < x_range[1]) - 1
        if not isinstance(states, list): states = [states]
        fig = plt.figure()
        ax = plt.axes()
        for state in states:
            ax.plot(self.results['grid'][x_min:x_max],
                    (self.results['wfns'][:,state][x_min:x_max]*scale)+self.results['energies'][state],
                    label=r'$\psi_{%s}$' %state)
        if on_pot and len(pot)>0:
            if self.in_AU:
                ax.plot(self.results['grid'][x_min:x_max], pot[x_min:x_max], label="potential")
            else:
                ax.plot(self.results['grid'][x_min:x_max],
                        ut.Constants.convert(pot[x_min:x_max], self.energy_units, to_AU=False),
                        label="potential")
        ax.set_xlim(x_range)
        if len(y_range) > 0:
            ax.set_ylim(y_range)
        #ax.set_xlabel(r'$\Delta r$')
        plt.legend()
        if save_file:
            plt.savefig(save_file)
        plt.show()


    def get_zpe(self):
        return self.results['energies'][0]





if __name__ == "__main__":
    import energy_arrays
    import utilities as ut
    import matplotlib.pyplot as plt

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
    gaussian_results = np.array([gaussian_grid, gaussian_energies])


    pot_mat = energy_arrays.ScanPotMat(grid=grid, gaussian_results=gaussian_results, min_shift=True)
    potential_matrix = pot_mat.matrix"""
    CH3_mass = ut.Constants.reduced_mass('C-H-H-H', to_AU=True)
    grid = Grid((0, 2 * np.pi), 101, inclusive=False).grid

    kin_mat = energy_arrays.CM_1D_kin_mat(grid, interval='0_to_2pi', mass=0.5)
    kinetic_matrix = kin_mat.matrix
    #plt.matshow(kinetic_matrix)

    pot_opts = {'type': 'periodic',
                'alpha':ut.Constants.convert(200, 'wavenumbers', to_AU=True)}
    pot_mat = energy_arrays.PredefinedPotMats(grid=grid, pot_opts=pot_opts)
    #plt.plot(grid, pot_mat.array)
    potential_matrix = pot_mat.matrix
    #potential_matrix = np.diag(np.zeros(len(grid)))

    file = '/Users/coire/McCoy/QOOH_repo/coire/DVR_dev/water_monomer_rpath/DVR_test'
    test = DVR_1D(grid=grid, kinetic_matrix=kinetic_matrix, potential_matrix=potential_matrix, file_name=file)
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
    test = DVR_1D(grid=grid, kinetic_matrix=kinetic_matrix, potential_matrix=potential_matrix, file_name=file)
    test.run()
    print("hello")



