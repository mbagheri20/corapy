import phonopy
import os
import subprocess
from pymatgen.core import Structure
import numpy as np
import re
import shutil
from spglib import get_symmetry, get_pointgroup, get_symmetry_dataset
from phonopy.structure.cells import estimate_supercell_matrix
from phonopy.interface.calculator import read_crystal_structure
from read_config import configs_reader
from pymatgen.io.vasp.outputs import Vasprun


configs = configs_reader('config.yaml')
ncore = configs[1]['cores']

current_dir = os.getcwd()
bs = Vasprun(current_dir+"/relaxation/vasprun.xml").as_dict()
egap = bs['output']['bandgap']

os.mkdir(current_dir+"/forces")
os.chdir(current_dir+"/forces")

shutil.copy(current_dir+"/relaxation/INCAR", current_dir+"/forces")
shutil.copy(current_dir+"/relaxation/POTCAR", current_dir+"/forces")
shutil.copy(current_dir+"/relaxation/CONTCAR", current_dir+"/forces")
shutil.copy(current_dir+"/relaxation/job_relax.sh", current_dir+"/forces")
os.rename('CONTCAR','POSCAR')
shutil.copy(current_dir+"/relaxation/CONTCAR", current_dir+"/forces")
os.rename('CONTCAR','POSCAR-unitcell')

struct = Structure.from_file('POSCAR')

# Find sc matrix
uc, _ = read_crystal_structure("POSCAR", interface_mode='vasp')
uc_cell = (uc.cell, uc.scaled_positions, uc.numbers)
dataset = get_symmetry_dataset(uc_cell, symprec=1e-5, angle_tolerance=-1.0, hall_number=0)
nmk = estimate_supercell_matrix(dataset,max_num_atoms=150)


# 1.2 generate supercell file
dim = "\""+str(nmk[0])+" "+ str(nmk[1])+ " "+ str(nmk[2])+ "\""
subprocess.run('phonopy -d --dim='+dim, shell=True)

# 2. Calculation of sets of forces
# 2.2 Modify INCAR
with open("POTCAR","r") as f:
    text = f.readlines()

ENs = []
for line in text:
    result = line.startswith("   ENMAX")
    if result is True:
        stripped = line.split(';', 1)[0]
        g = re.findall("\d+\.\d+", stripped)
        ENs.append((float(g[0])))

ENMAX = max(ENs, key=lambda x: float(x))


myincar ="PREC = Accurate\n"\
         +"GGA = PS\n"\
         +"IBRION = -1\n"\
         +"NELMIN = 5\n"\
         +"ENCUT = "+str(format(1.3*ENMAX, ".6f"))+"\n"\
         +"EDIFF = 1.000000e-08\n"\
         +"ISMEAR = 0\n"\
         +"SIGMA = 1.000000e-02\n"\
         +"IALGO = 38\n"\
         +"LREAL = .FALSE.\n"\
         +"ADDGRID = .TRUE.\n"\
         +"LWAVE = .FALSE.\n"\
         +"LCHARG = .FALSE.\n"\
         +"NPAR = " + str(int(ncore/10))

#open('INCAR', 'wt').write(myincar)
with open("INCAR","r") as f:
    incar_relax = f.read()

incar_relax = incar_relax.replace('NSW = 100\n','')
incar_relax = incar_relax.replace('ISIF = 3\n','')
incar_relax = incar_relax.replace('EDIFFG = -1.000000e-06\n','')
incar_force = incar_relax.replace('IBRION = 2\n','IBRION = -1\n')

open('INCAR', 'wt').write(incar_force)

struct_sc = Structure.from_file('SPOSCAR')
if egap > 2:
    Rk=20
elif egap < 1:
    Rk=40
else:
    Rk=30


N1 = int(max(1,Rk*struct_sc.lattice.reciprocal_lattice.a/(2*np.pi)+0.5))
N2 = int(max(1,Rk*struct_sc.lattice.reciprocal_lattice.b/(2*np.pi)+0.5))
N3 = int(max(1,Rk*struct_sc.lattice.reciprocal_lattice.c/(2*np.pi)+0.5))
mykpt = 'Automatic mesh\n'+'0\n'+'Gamma\n'+ str(N1)+" "+str(N2)+" "+str(N3)+'\n'+'0 0 0'
fout = open("KPOINTS", 'wt').write(mykpt)

os.rename('job_relax.sh','job_force.sh')



path = current_dir+"/forces"
npos=0
for file in os.listdir(path):
        if file.startswith("POSCAR-0"):
            npos=npos+1

foutt = open("npos.txt", "wt").write(str(npos))

# Displacements
for i in range(1, npos+1):
    dir_name = f"disp-{i:03d}"
    os.mkdir(dir_name)
    files_to_copy = ["KPOINTS", "POTCAR", "INCAR", "job_force.sh"]
    for file_name in files_to_copy:
        shutil.copy(file_name, dir_name)

    # Rename and move POSCAR file
    poscar_src = f"POSCAR-{i:03d}"
    poscar_dest = os.path.join(dir_name, "POSCAR")
    os.rename(poscar_src, poscar_dest)

for i in range(1, npos+1):
    dir_name = f"disp-{i:03d}"
    os.chdir(dir_name)
    subprocess.run("sbatch job_force.sh", shell=True)

    os.chdir("..")



os.remove('job_force.sh')

