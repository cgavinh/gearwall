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
                 num_points):
        self.domain = domain
        self.num_points = num_points
        self.initialize()

    def initialize(self):
        self.grid = np.linspace(self.domain[0], self.domain[1], self.num_points)

class AnalyzeDVR:
    def __init__(self,
                 results_file,
                 file_type='pickle',
                 in_AU = True,
                 energy_units='wavenumbers',
                 length_units='angstroms'):
        self.results_file = results_file
        self.file_type = file_type
        self._results = None
        self.in_AU = in_AU
        self.energy_units=energy_units
        self.length_units=length_units
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
        self.results['grid'] = ut.Constants.convert(self.results['grid'], self.length_units, to_AU=False)
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

    def plot_wfns(self, states, pot = [], on_pot=False, range=[]):
        if not len(range) > 1:
            range = [self.results['grid'].min(), self.results['grid'].max()]
        if not isinstance(states, list): states = [states]
        fig = plt.figure()
        ax = plt.axes()
        for state in states:
            ax.plot(self.results['grid'], self.results['wfns'][:,state]+self.results['energies'][state], label=r'$\psi_{%s}$' %state)
        if on_pot and len(pot)>0:
            if self.in_AU:
                ax.plot(self.results['grid'], pot, label="potential")
            else:
                ax.plot(self.results['grid'],
                        ut.Constants.convert(pot, self.energy_units, to_AU=False),
                        label="potential")
        ax.set_xlim(range)
        ax.set_xlabel(r'$\Delta r$')
        plt.legend()
        plt.show()

    def get_zpe(self, in_AU = True, energy_units=''):
        return self.results['energies'][0]





if __name__ == "__main__":
    import array_tools
    import utilities as ut
    import matplotlib.pyplot as plt

    OH_mass = ut.Constants.reduced_mass('O-H', to_AU=True)
    grid_range = ut.Constants.convert(1, 'angstroms', to_AU=True)
    freq = ut.Constants.convert(3300, 'wavenumbers', to_AU=True)

    grid = Grid((-grid_range,grid_range), 100).grid
    kin_mat = array_tools.CM_1D_kin_mat(grid, interval='-infty_to_infty', mass=OH_mass)
    kinetic_matrix = kin_mat.matrix

    pot_opts = {'type': 'harmonic oscillator',
                'mass': OH_mass,
                'frequency': freq}

    pot_mat = array_tools.AnalyticalPotMats(grid, pot_opts)
    potential_matrix = pot_mat.matrix

    file = '/Users/coire/McCoy/QOOH_repo/coire/DVR_dev/DVR_test'
    test = DVR_1D(grid=grid, kinetic_matrix=kinetic_matrix, potential_matrix=potential_matrix, file_name=file)
    test.run()
    anl_test = AnalyzeDVR(results_file=file, in_AU=False, energy_units='wavenumbers', length_units='angstroms')
    anl_test.get_frequency((0,1))
    anl_test.rephase_wfns([0])
    anl_test.plot_wfns(states=[0,1,2],
                       pot=np.diag(potential_matrix),
                       on_pot=True)
    zpe = anl_test.get_zpe()
    print("hello")



