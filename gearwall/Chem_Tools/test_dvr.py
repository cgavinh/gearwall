# A space to test DVR code
from Chem_Tools import dvr, gaussian_tools as gt, chem_utilities as uts, dvr_arrays
import numpy as np
from scipy import interpolate

def test_infty_to_infty_1d_DVR(range_angstroms=(0.5, 1.5), re_angstroms=1, atom_str='O-H', num_grid_pts=50, freq_cm=3300):
    range_bohr = ([uts.Constants.convert(x, 'angstroms', to_AU=True) for x in range_angstroms])
    freq_ha = uts.Constants.convert(freq_cm, 'wavenumbers', to_AU=True)
    grid_bohr = dvr.Grid(domain=range_bohr, num_points=num_grid_pts, inclusive=True).grid
    mass_me = uts.Constants.reduced_mass(atom_str, to_AU=True)
    re_bohr = uts.Constants.convert(re_angstroms, 'angstroms', to_AU=True)

    pot_opts = {'type': 'harmonic oscillator',
                'mass': mass_me,
                'frequency': freq_ha,
                're': re_bohr}
    potential_matrix = dvr_arrays.PredefinedPotMats(grid=grid_bohr,
                                                    pot_opts=pot_opts)
    kinetic_matrix = dvr_arrays.CM_1D_kin_mat(grid=grid_bohr,
                                              interval='-infty_to_infty',
                                              mass=mass_me)
    dvr_test = dvr.DVR_1D(grid=grid_bohr,
                          kinetic_matrix=kinetic_matrix.matrix,
                          potential_matrix=potential_matrix.matrix,
                          filename="test_infty_to_infty_1d_DVR",
                          save_method='pickle')
    dvr_test.run()
    dvr_test_a = dvr.AnalyzeDVR(results_file="test_infty_to_infty_1d_DVR",
                                save_method='pickle',
                                in_AU=True,
                                energy_units='wavenumbers',
                                grid_unit='angstroms')
    dvr_test_a.plot_wfns(states=[0,1,2,3,4],
                         on_pot=True,
                         save_file='test_infty_to_infty_1d_DVR_wfns.png',
                         scale=20000)
    dvr_test_a.plot_energies(energies_range=[0,20], energy_step=2)
    dvr_test_a.plot_ind_wfns(scale=30000)
    a=4

if __name__ == "__main__":
    #test_infty_to_infty_1d_DVR()
    #DVR()
    #print("hello")
    clhocl_oop_s = gt.GLogInterpreter("../hocl/oop_s_rOH_opt/clhocl/clhocl_oop_s.log")
    grid = dvr.Grid(domain=[-75 * 2 * np.pi / 360, 75 * 2 * np.pi / 360], num_points=100).grid

    f_linear = interpolate.interp1d(np.linspace(-75 * 2 * np.pi / 360, 75 * 2 * np.pi / 360, 76),
                                    clhocl_oop_s.Energy,
                                    kind='linear', fill_value='extrapolate',
                                    bounds_error=False)  # linear spline

    pot_array = f_linear(grid)
    pot_array = pot_array - pot_array.min()
    pot_matrix = np.diag(pot_array)

    # create the mass grid
    ## first I need to evaluate rOH at the grid points
    f_rOH = interpolate.interp1d(np.linspace(-75 * 2 * np.pi / 360, 75 * 2 * np.pi / 360, 76),
                                 uts.Constants.convert(clhocl_oop_s.intCoords['r3'], to_AU=True, unit="angstroms"),
                                 kind='linear', fill_value='extrapolate',
                                 bounds_error=False)
    rOH_grid = f_rOH(grid)
    ## now I need to construct the mass matrix
    mass_matrix = np.zeros((len(grid), len(grid)))


    def g(x):
        return uts.Constants.reduced_mass("O-H") * rOH_grid[x] ** 2


    for i in range(len(mass_matrix)):
        for j in range(len(mass_matrix)):
            mass_matrix[i][j] = 0.5 * (g(i) + g(j))

    kinetic = dvr_arrays.CM_1D_kin_mat(grid=grid,
                                       interval='-infty_to_infty',
                                       mass=mass_matrix)

    clhocl_oop_dvr = dvr.DVR_1D(grid=grid, kinetic_matrix=kinetic.matrix, potential_matrix=pot_matrix,
                                filename="clhocl_oop_dvr")

    clhocl_oop_dvr.run()

    clhocl_oop_anDVR = dvr.AnalyzeDVR("clhocl_oop_dvr", grid_unit='radians')

    clhocl_oop_anDVR.plot_wfns(states=[0, 1, 2], on_pot=True, scale=10000, save_file="clhocl_oop_dvr_wfns.pdf")
    print("hello")