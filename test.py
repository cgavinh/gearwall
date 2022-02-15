import density_cube
import gaussian_tools as gt
import glob
from string import Template
import sys
import misc
import os
sys.path.append('../mccoy_hoono/')
import numpy as np
import McUtils.GaussianInterface as GI

if __name__ == "__main__":
    atoms = ['H', 'O', 'O', 'N', 'O']
    specs = np.array([[1, 2, 3, 4], [1, 2, 3], [1, 2]])
    names = np.array([['r1', 'r2', 'r3', 'r4'], ['a1', 'a2', 'a3'], ['d1', 'd2']])
    vals = np.array([[0.96897556, 1.45491931, 1.45878406, 1.17072885],
                     [97.59494903, 102.29752444, 109.06955123],
                     [180., 180.]])
    units = ["Angstroms", "Degrees"]
    hoono_eq = misc.zmat(atoms=atoms, easy_build_dict={'specs':specs, 'names':names, 'vals':vals, 'units':units})
    # hoono_eq_log = gt.GLogInterpreter(log_files='../mccoy_hoono/hoono_eq.log')
    # hoono_eq_log.pull_forces()
    # params = {
    #                 'header': f'%nproc=28\n%mem=120GB\n%chk=',
    #                 'job': '#p SP mp2/aug-cc-pvtz density=current\n\n',
    #                 'description': 'title\n\n',
    #                 'charge': '0 1\n',
    #                 'molecule': '',
    #                 'variables': 'Variables\n',
    #                 'vars': '',
    #                 'whitespace': '\n'*10}
    # disps = misc.target_displacements(eq_zmat=hoono_eq,
    #                           force_constants=hoono_eq_log.force_constants)
    # misc.write_displaced_zmats(eq_zmat=hoono_eq,
    #                            force_constants=hoono_eq_log.force_constants,
    #                            fname_base='hoono',
    #                            params=params,
    #                            job_file_path='../mccoy_hoono/explore_disps/')

    # hoono_fchk = gt.FchkInterpreter(fchks=['../mccoy_hoono/explore_disps/hoono_a1_95p93353.fchk'])
    # hoono_log = gt.GLogInterpreter(log_files=['../mccoy_hoono/explore_disps/hoono_a1_95p93353.log'])
    # with GI.GaussianLogReader('../mccoy_hoono/hoono_eq.log') as reader:
    #     parse = reader.parse("OptimizedScanEnergies")
    # z = misc.zmat(g_str='\n Charge =  0 Multiplicity = 1\n H\n O                    1    r1\n O                    2    r2       1    a1\n N                    3    r3       2    a2       1    d1       0\n O                    4    r4       3    a3       2    d2       0\n       Variables:\n  r1                    0.96898                  \n  r2                    1.45492                  \n  r3                    1.45878                  \n  r4                    1.17073                  \n  a1                   95.93353                  \n  a2                  102.29752                  \n  a3                  109.06955                  \n  d1                  180.                       \n  d2                  180.')
    #misc.zmat()
    # result = gt.GaussianResults(files=['../mccoy_hoono/explore_disps/hoono_d2_176p1845'])
    # result = gt.GaussianResults(files=[os.path.splitext(val)[0] for val in glob.glob("../mccoy_hoono/explore_disps/*.gjf")])
    # files=['../mccoy_hoono/explore_disps/hoono_r3_1p437917',
    #        '../mccoy_hoono/explore_disps/hoono_a2_100p8428',
    #        '../mccoy_hoono/explore_disps/hoono_d1_188p7757',
    #        '../mccoy_hoono/explore_disps/hoono_a3_110p1829',
    #        '../mccoy_hoono/explore_disps/hoono_r4_1p162801',
    #        '../mccoy_hoono/explore_disps/hoono_r4_1p178656',
    #        '../mccoy_hoono/explore_disps/hoono_a1_99p25636',
    #        '../mccoy_hoono/explore_disps/hoono_d1_171p2242',
    #        '../mccoy_hoono/explore_disps/hoono_a2_103p7522',
    #        '../mccoy_hoono/explore_disps/hoono_r1_0p958995',
    #        '../mccoy_hoono/explore_disps/hoono_r1_0p978956',
    #        '../mccoy_hoono/explore_disps/hoono_a3_107p9562',
    #        '../mccoy_hoono/explore_disps/hoono_r3_1p479650',
    #        '../mccoy_hoono/explore_disps/hoono_a1_95p93353',
    #        '../mccoy_hoono/explore_disps/hoono_r2_1p441450',
    #        '../mccoy_hoono/explore_disps/hoono_r2_1p468388']
    # hoono_eq_energy = -280.4467347
    # results = misc.pull_displaced_results(eq_zmat=hoono_eq,
    #                                       eq_energy=hoono_eq_energy,
    #                                       files=files,
    #                                       sec_der=True)
    # result = gt.GaussianResults(files=[os.path.splitext(val)[0] for val in glob.glob("../mccoy_hoono/10cm_disp/*.gjf")])
    # result = gt.GaussianResults(files=[os.path.splitext(val)[0] for val in glob.glob("../mccoy_hoono/30cm_disp/*.gjf")])
    with GI.GaussianLogReader('../mccoy_hoono/explore_disps/hoono_d2_176p1845.log') as reader:
        parse = reader.parse("InputZMatrix")
    print("Hello")
