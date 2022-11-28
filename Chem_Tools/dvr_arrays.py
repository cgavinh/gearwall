import numpy as np
from Chem_Tools import gaussian_tools as gt


class CM_1D_kin_mat:
    def __init__(self,
                 grid,
                 interval,
                 mass=1):
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
                    return ((-1) ** (i - j)) * (np.power(np.pi, 2) / 3)
                else:
                    return ((-1) ** (i - j)) * (2 / np.power((i - j), 2))
            self.matrix =  coeff * self.fill_array(kinE, len(self.grid), len(self.grid))

        if self.interval =='0_to_2pi':
            N = (len(self.grid)-1)/2
            coeff = 1/(2.*self.mass)
            for i in range(len(self.matrix)):
                for j in range(len(self.matrix)):
                    if i == j:
                        self.matrix[i][j] = coeff*(N*(N+1.))/3.
                    else:
                        self.matrix[i][j] = coeff*((-1)**(i -j))*(np.cos((np.pi*(i-j))/(2.*N+1.))/(2.*(np.sin((np.pi*(i-j))/(2.*N+1.)))**2.))


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

class PredefinedPotMats:
    """
    for HO, pot_opts should include (type='harmonic oscillator', mass, frequency, and re)
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
            self.array, self.matrix = self.HO()
        if self.pot_opts['type'] == 'periodic':
            self.array, self.matrix = self.periodic()

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
            pot_grid[i] = (1 / 2) * self.pot_opts['mass'] * self.pot_opts['frequency'] ** 2 * (self.grid[i]-self.pot_opts['re']) ** 2
        pot_matrix = np.diag(pot_grid)
        return pot_grid, pot_matrix

    def periodic(self):
        pot_grid = np.zeros(self.grid.shape)
        for i in range(len(self.grid)):
            pot_grid[i] = self.pot_opts['alpha']*(1+np.cos(3*self.grid[i]))
        pot_matrix = np.diag(pot_grid)
        return pot_grid, pot_matrix

class GaussianPots:
    """"""
    def __init__(self,
                 files,
                 grid):
        self.grid = grid
        self.frames = gt.GLogInterpreter(files)

"""class ScanPotMat:
    def __init__(self,
                 grid,
                 scan_results,
                 min_shift=True,
                 int_opts):
        self.grid = grid
        self.scan_results = scan_results
        self.min_shift = min_shift
        self.int_opts = int_opts
        self.initialize()


    def initialize(self):
        f_cubic = interpolate.interp1d(self.scan_results[0],
                                       self.scan_results[1],
                                       kind='cubic', fill_value=np.nan,
                                       bounds_error=False)  # cubic spline
        f_linear = interpolate.interp1d(self.scan_results[0],
                                        self.scan_results[1],
                                        kind='linear', fill_value='extrapolate',
                                        bounds_error=False)  # linear spline
        if self.int_opts['ext']:
            if self.int_opts['in_range'] == 'cubic':


            self.pot_array  = np.where(np.isnan(f_cubic(self.grid)), f_linear(self.grid), f_cubic(self.grid))
        if self.min_shift:
            self.pot_array = self.pot_array - np.min(self.pot_array)
        self.matrix =  np.diag(self.pot_array)

    def plot_pot(self):
        fig = plt.figure()
        axes = plt.axes()
        plt.scatter(self.scan_results[0], self.scan_results[1], label="input data")
        plt.plot(self.grid, self.pot_array, label="interpolated potential")
        plt.legend()
        plt.show()"""
