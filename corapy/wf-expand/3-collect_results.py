import phonopy
import os
import requests
from ase.io import read
from ase import io,Atoms
import subprocess
from pymatgen.core import Structure
import numpy as np
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.sets import MPRelaxSet
import re
import shutil
import sys
from spglib import get_symmetry, get_pointgroup, get_symmetry_dataset
from phonopy.structure.cells import estimate_supercell_matrix
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure


current_dir = os.getcwd()
os.mkdir(current_dir+"/results")
os.chdir(current_dir+"/forces")

npos =open(current_dir+"/forces/""npos.txt", "rt").read()
phc = 'phonopy -f disp-{001..'+npos+'}/vasprun.xml'
path_to_bash = "/bin/bash"  # or whatever is appropriate
process = subprocess.Popen(phc,
                           stdout=subprocess.PIPE,
                           shell=True,
                           executable=path_to_bash)
output, error = process.communicate()



shutil.copy(current_dir+"/forces/FORCE_SETS",current_dir+"/results")
shutil.copy(current_dir+"/forces/POSCAR-unitcell",current_dir+"/results")
shutil.copy(current_dir+"/forces/phonopy_disp.yaml",current_dir+"/results/disp.yaml")
shutil.copy(current_dir+"/relaxation/INCAR",current_dir+"/results/INCAR-relax")
shutil.copy(current_dir+"/relaxation/KPOINTS",current_dir+"/results/KPOINTS-relax")
shutil.copy(current_dir+"/forces/disp-001/INCAR",current_dir+"/results/INCAR-force")
shutil.copy(current_dir+"/forces/disp-001/KPOINTS",current_dir+"/results/KPOINTS-force")


# make phonon.yaml file


#uc = read('POSCAR')
uc, _ = read_crystal_structure("POSCAR", interface_mode='vasp')
dataset = get_symmetry_dataset(uc, symprec=1e-5, angle_tolerance=-1.0, hall_number=0)
nmk = estimate_supercell_matrix(dataset,max_num_atoms=150)
phonon = Phonopy(uc,[[nmk[0], 0, 0], [0, nmk[1], 0], [0, 0, nmk[2]]],primitive_matrix='auto')
phonon.save(filename="phonon.yaml")
shutil.copy(current_dir+"/forces/phonon.yaml",current_dir+"/results")
