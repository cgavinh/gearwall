import sys
sys.path.append("/Users/coire/McCoy/") #Path to McUtils, References..etc
import numpy as np
import copy
from McUtils import Data as McData
from Chem_Tools import gaussian_tools as gt


def fourth_ders(results):
    der = {'disp_vals': [],
           'fourth_ders': []}

    derivatives = {}

    for q in results:
        new_q = copy.deepcopy(der)
        eq_index = int((len(results[q]['disp_vals']) - 1) / 2)
        for i in range(1, eq_index + 1):
            h = results[q]['disp_vals'][eq_index + i]
            new_q['disp_vals'].append(h)
            new_q['fourth_ders'].append((results[q]['force_constants'][eq_index + i] +results[q]['force_constants'][eq_index - i] -2 * results[q]['force_constants'][eq_index]) /h**2)
        derivatives[q] = new_q
    return derivatives


def pull_displaced_results(eq_zmat, eq_energy, files, eq_log):
    '''
    eq_energy in hartree
    outputs displaced energy in wavenumbers'''
    r = {'filenames': [],
         'disp_vals': [],
         'disp_units': [],
         'disp_energies': [],
         'force_constants': []
         }
    results = {}
    gr = gt.GaussianResults(files = files)
    gr.pull_forces()
    #results['filenames'] = files
    #results['disp_energies'] = gr.MP2Energy
    for i in range(len(files)):
        #find displaced coord
        for e in range(len(eq_zmat.coords)):
            if round(gr.input_zmat[i].coords[e].val,5) != round(eq_zmat.coords[e].val,5):
                if gr.input_zmat[i].coords[e].name in results:
                    results[gr.input_zmat[i].coords[e].name]['filenames'].append(files[i])
                    results[gr.input_zmat[i].coords[e].name]['disp_vals'].append(gr.input_zmat[i].coords[e].val - eq_zmat.coords[e].val)
                    results[gr.input_zmat[i].coords[e].name]['disp_units'].append(gr.input_zmat[i].coords[e].units)
                    results[gr.input_zmat[i].coords[e].name]['disp_energies'].append((gr.MP2Energy[i]-eq_energy)*McData.UnitsData.convert( "Hartrees","Wavenumbers"))
                    results[gr.input_zmat[i].coords[e].name]['force_constants'].append(gr.force_constants[i][e][e])
                else:
                    new_result = copy.deepcopy(r)
                    new_result['filenames'].append(files[i])
                    new_result['disp_vals'].append(gr.input_zmat[i].coords[e].val - eq_zmat.coords[e].val)
                    new_result['disp_units'].append(gr.input_zmat[i].coords[e].units)
                    new_result['disp_energies'].append((gr.MP2Energy[i]-eq_energy)*McData.UnitsData.convert( "Hartrees","Wavenumbers"))
                    new_result['force_constants'].append(gr.force_constants[i][e][e])
                    results[gr.input_zmat[i].coords[e].name] = new_result
    #Add equilibrium geometry and convert force units
    eq = gt.GLogInterpreter(log_files=[eq_log + '.log'])
    eq.pull_forces()
    forces = convert_force_units(eq_zmat=eq_zmat, force_constants=eq.force_constants)
    for r in results:
        results[r]['filenames'].append("eq")
        results[r]['disp_vals'].append(0)
        results[r]['disp_units'].append(results[r]['disp_units'][-1])
        results[r]['disp_energies'].append(0)
        results[r]['force_constants'] = [a * McData.UnitsData.convert("Hartrees", "Wavenumbers") for a in
                                         results[r]['force_constants']]
        if results[r]['disp_units'][0] ==  "Angstroms":
            results[r]['force_constants'] = [a * (1/McData.UnitsData.convert("BohrRadius", "Angstroms"))**2 for a in results[r]['force_constants']]
        elif results[r]['disp_units'][0] ==  "Degrees":
            results[r]['force_constants'] = [a * ((2 * np.pi) / 360)**2 for a in results[r]['force_constants']]
        for e in range(len(eq_zmat.coords)):
            if eq_zmat.coords[e].name == r:
                results[r]['force_constants'].append(forces[e][e])

    for key in results:
        zipped_lists = zip(results[key]['disp_vals'],
                           results[key]['filenames'],
                           results[key]['disp_units'],
                           results[key]['disp_energies'],
                           results[key]['force_constants'])
        sorted_pairs = sorted(zipped_lists) #sorts by items in the first list
        tuples = zip(*sorted_pairs)
        results[key]['disp_vals'], results[key]['filenames'],results[key]['disp_units'],results[key]['disp_energies'],results[key]['force_constants'] = [list(tuple) for tuple in tuples]

    # if sec_der:
    #     for key in results:
    #         results[key]['2nd_der'] = (results[key]['disp_energies'][0] + results[key]['disp_energies'][1])/results[key]['disp_vals'][0]**2
    return results



def convert_force_units(eq_zmat, force_constants):
    force_constants = force_constants * McData.UnitsData.convert("Hartrees", "Wavenumbers")
    for i in range(len(force_constants)):
        if eq_zmat.coords[i].units == "Angstroms":
            force_constants[i] = force_constants[i] * (1 / McData.UnitsData.convert("BohrRadius", "Angstroms"))
            force_constants[:, i] = force_constants[:, i] * (1 / McData.UnitsData.convert("BohrRadius", "Angstroms"))
        elif eq_zmat.coords[i].units == "Degrees":
            force_constants[i] = force_constants[i] * (2 * np.pi) / 360
            force_constants[:, i] = force_constants[:, i] * (2 * np.pi) / 360
    return force_constants

def target_displacements(eq_zmat, force_constants, disp_en = 20):
    #Gaussian force constants are in Hartree/Bohr^2 or Hartree/Radian^2
    force_constants = convert_force_units(eq_zmat = eq_zmat, force_constants=force_constants)
    disps = np.zeros(len(eq_zmat.coords))
    for d in range(len(disps)):
        disps[d] = np.sqrt(abs(2 * disp_en / force_constants[d][d]))
    return disps

def write_displaced_zmats(eq_zmat,
                          force_constants,
                          disp_en = 20,
                          sbatch_params='',
                          fname_base='',
                          params={},
                          job_file_path='',
                          fchks=False):
    with open(job_file_path + f'g16_{fname_base}_{disp_en}.sh', 'w', newline='\n') as sbatch:
        sbatch.write('#!/bin/bash\n')
        sbatch.write(sbatch_params)
        sbatch.write('\n# load Gaussian environment\nmodule load contrib/g16.b01\n\n# '
                     'debugging information\necho '
                     '\"**** Job Debugging Information ****\"\necho \"This job will '
                     'run on $SLURM_JOB_NODELIST"\necho \"\"\n'
                     'echo "ENVIRONMENT VARIABLES"\nset\necho \"********************'
                     '**************************\"\n\n# run Gaussian\n')
        disps = target_displacements(eq_zmat=eq_zmat, force_constants=force_constants, disp_en=disp_en)
        for d in range(len(disps)):
            pm = [disps[d], -disps[d]]
            for i in range(2):
                new_zmat = copy.deepcopy(eq_zmat)
                new_zmat.coords[d].val += pm[i]
                if len(str(new_zmat.coords[d].val)) > 8:
                    fname = fname_base + f"_{new_zmat.coords[d].name}_{str(new_zmat.coords[d].val).replace('.', 'p')[0:8]}"
                else:
                    fname = fname_base + f"_{new_zmat.coords[d].name}_{str(new_zmat.coords[d].val).replace('.', 'p')}"
                job_params = copy.deepcopy(params)
                job_params['header'] += f'{fname}.chk\n'
                job_params['molecule'] = new_zmat.print_mat()
                job_params['vars'] = new_zmat.print_vals()
                job = gt.GaussianJob(params=job_params, filename=f'{job_file_path}{fname}.gjf')
                job.save_job()
                sbatch.write(f'g16 {fname}.gjf\n')
                if fchks:
                    sbatch.write(f'formchk {fname}.chk\n')




class int_coord:
    def __init__(self, name=None, spec=None, val=None, units=None):
        self.name=name
        self.spec=spec
        self.val=val
        self.units=units


class zmat:

    #Example Build Dict
    # build_dict={"atoms" : ['H', 'O', 'N', 'O'],
    #             "specs" :  np.array([[1, 2, 3], [1, 2], [1]]),
    #             "names" : np.array([['r1', 'r2', 'r3'], ['a1', 'a2'], ['d3']]),
    #             "vals" : np.array([[0.97026274, 1.42776435, 1.17680217],
    #                                [101.62778541, 110.70859147],
    #                                [180.]]),
    #             "units" : ["Angstroms", "Degrees"]}

    def __init__(self, atoms=None, coords=None, charge=[0,1], easy_build_dict=None, g_str=None):
        self.atoms=atoms
        self.coords=coords
        self.charge=charge
        self.build_dict=easy_build_dict
        self.g_str = g_str
        if self.coords is None and self.build_dict is not None:
            self.coords=[]
            specs = []# build #N-6 array that specifies all connections
            for i in range(len(self.build_dict['specs'])):
                for j in range(len(self.build_dict['specs'][i])):
                    new_spec = []
                    new_spec.append(j+i+2)
                    for k in range(i+1):
                        new_spec.append(self.build_dict['specs'][k][j+i-k])
                    specs.append(new_spec)
            i=0
            u = 0
            for a in range(len(self.build_dict['names'])):
                for b in range(len(self.build_dict['names'][a])):
                    if len(specs[i]) > 2:
                        u = 1
                    else:
                        u = 0
                    self.coords.append(int_coord(name=self.build_dict['names'][a][b],
                                                 spec=specs[i],
                                                 val=self.build_dict['vals'][a][b],
                                                 units=self.build_dict['units'][u]))
                    i += 1
        elif self.coords is None and self.g_str is not None:
            self.coords = []
            self.atoms = []
            lines = g_str.splitlines()
            i = 0
            n = 0
            b = []
            a = []
            d = []

            for line in lines:
                if 'Charge' in line:
                    self.charge[0] = line.split()[2]
                    self.charge[1] = line.split()[5]
                    i = 1
                    continue
                elif 'Variables' in line:
                    i = 2
                    for bond in b:
                        self.coords.append(bond)
                    for ang in a:
                        self.coords.append(ang)
                    for dih in d:
                        self.coords.append(dih)
                    continue
                elif i == 1:
                    n += 1
                    l = line.split()
                    self.atoms.append(l[0])
                    if len(l) == 3:
                        b.append(int_coord(name=l[2], spec = [n, int(l[1])], units="Angstroms"))
                    elif len(l) == 5:
                        b.append(int_coord(name=l[2], spec=[n, int(l[1])], units="Angstroms"))
                        a.append(int_coord(name=l[4], spec=[n, int(l[1]), int(l[3])], units="Degrees"))
                    elif len(l) > 5:
                        b.append(int_coord(name=l[2], spec=[n, int(l[1])], units="Angstroms"))
                        a.append(int_coord(name=l[4], spec=[n, int(l[1]), int(l[3])], units="Degrees"))
                        d.append(int_coord(name=l[6], spec=[n, int(l[1]), int(l[3]), int(l[5])], units="Degrees"))
                elif i == 2:
                    l = line.split()
                    for c in self.coords:
                        if c.name == l[0]:
                            c.val = float(l[1])
        elif self.coords is None and self.g_str is None:
            print(f'self.coords is None and self.g_str is None')
            print(self.atoms)






    def print_mat(self):
        # z_mat = ''
        # z_mat += self.atoms[0] + '\n'
        # if len(self.atoms) > 1:
        #     z_mat += f"{self.atoms[1]} {self.coords[0].spec[-1]} {self.coords[0].name}" + '\n'
        # if len(self.atoms) > 2:
        #     z_mat += f"{self.atoms[2]} {self.coords[1].spec[-1]} {self.coords[1].name} {self.refs[2][0]} {self.vars[1][0]}" + '\n'
        # if len(self.atoms) > 3:
        #     for i in range(3, len(self.refs[0])):
        #         z_mat += f"{self.refs[0][i]} {self.refs[1][i - 1]} {self.vars[0][i - 1]} {self.refs[2][i - 2]} {self.vars[1][i - 2]} {self.refs[3][i - 3]} {self.vars[2][i - 3]}" + '\n'
        # return z_mat

        #assume order of coords is all bonds, all angles, all dihedrals as they appear on z-mat
        pre_mat = [[''],[''],['',''],['',''],['','',''],['','','']]
        for coord in self.coords:
            if len(coord.spec) == 2:
                pre_mat[0].append(coord.spec[-1])
                pre_mat[1].append(coord.name)
            elif len(coord.spec) == 3:
                pre_mat[2].append(coord.spec[-1])
                pre_mat[3].append(coord.name)
            elif len(coord.spec) == 4:
                pre_mat[4].append(coord.spec[-1])
                pre_mat[5].append(coord.name)

        mat = np.full((len(self.atoms),7), '', dtype='<U10') #weirdness with strings and data types. cuts off string if more than 10 chars long
        mat[:,0] = self.atoms
        for i in range(6):
            mat[:,i+1] = pre_mat[i]

        mat_str = ''
        for i in range(len(mat)):
            line=''
            for e in range(len(mat[i])):
                line += mat[i][e] + ' '
            mat_str += line + '\n'

        return mat_str

    def print_vals(self):
        vals = ''
        for coord in self.coords:
            vals += f' {coord.name}={coord.val}\n'
        return vals















