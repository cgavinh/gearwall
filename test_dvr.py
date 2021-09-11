# A space to test DVR code
import dvr
import dvr_arrays
import utilities as uts

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

#if __name__ == "__main__":
    #test_infty_to_infty_1d_DVR()
    #DVR()
    #print("hello")