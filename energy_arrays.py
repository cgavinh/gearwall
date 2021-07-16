import numpy as np
import utilities as ut

class CM_1D_kin_mat:
    def __init__(self,
                 grid,
                 interval,
                 mass):
        self.grid = grid
        self.interval=interval
        self.mass = mass
        self.intialize()

    def intialize(self):
        self.interval = self.interval.lower()
        self.matrix = np.zeros((len(self.grid), len(self.grid)))
        dx = self.grid[1] - self.grid[0]
        if self.interval =='-infty_to_infty':
            coeff = 1 / (2 * self.mass * dx ** 2)
            def kinE(i, j):
                if i == j:
                    return coeff * ((-1) ** (i - j)) * (np.power(np.pi, 2) / 3)
                else:
                    return coeff * ((-1) ** (i - j)) * (2 / np.power((i - j), 2))
            self.matrix =  self.fill_array(kinE, len(self.grid), len(self.grid))

        if self.interval =='0_to_infty':
            def kinE(i, j):
                if i == j:
                    return coeff * ((-1) ** (i - j)) * (np.power(np.pi, 2) / 3)
                else:
                    return coeff * ((-1) ** (i - j)) * (2 / np.power((i - j), 2))
            self.matrix = self.fill_array(kinE, len(self.grid), len(self.grid))

    def fill_array(self, expression, nrows, ncolumns):
        """
        :param expression: a function that takes two indices (i, j) and returns a value
        :param nrows: the number of rows
        :type nrows: int
        :param ncolumns: the number of columns
        :type ncolumns: int
        """
        my_array = np.zeros((nrows, ncolumns))
        for i in range(nrows):
            for j in range(ncolumns):
                my_array[i][j] = expression(i, j)
        return my_array

class AnalyticalPotMats:
    """
    for HO, pot_opts should include (type='harmonic oscillator', mass, and frequency)
    for MO, pot_opts should include (type='morse oscillator', alpha, dissociation_energy, re)
    """
    def __init__(self,
                 grid,
                 pot_opts):

        self.grid = grid
        self.pot_opts = pot_opts
        self.initialize()

    def initialize(self):
        if self.pot_opts['type'] == 'harmonic oscillator':
            self.matrix = self.HO()

    def HO(self):
        """Here we will use a harmonic oscillator function to create the matrix representation of the potential energy.
           :param grid: the grid created by dvr_grid (atomic units)
           :type grid: np.array
           :param mass: The reduced mass of your system (atomic units)
           *Note: this mass is the same as the mu you will use in the kinetic energy function, but since your potential will
           not always rely on mass so you do want to define them twice.
           :type mass: float
           :param frequency: the frequency of the vibration (atomic units)
           :type frequency: float
           :param re: Equilibrium distance of vibration (atomic units)
           :type re: float
           :return: Potential energy array
           :rtype: np.array """
        # First solve for the HO potential energy at each grid point
        # Then project your results into the matrix representation!
        # Remember: The potential matrix is 0 for all terms off the diagonal for any type of potential in DVR
        pot_grid = np.zeros(self.grid.shape)
        for i in range(len(self.grid)):
            pot_grid[i] = (1 / 2) * self.pot_opts['mass'] * self.pot_opts['frequency'] ** 2 * self.grid[i] ** 2
        pot_matrix = np.diag(pot_grid)
        return pot_matrix