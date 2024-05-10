#! /usr/bin/env python
import os
from pymatgen.core import Structure
import numpy as np
from mp_api.client import MPRester
from pymatgen.io.vasp.sets import MPRelaxSet
import re
import shutil
from pymatgen.io.vasp import Poscar
from read_config import configs_reader
from Outcar_checking import force_check
import subprocess
from pymatgen.io.vasp.outputs import Vasprun

# for legacy MP
#from pymatgen.ext.matproj import MPRester
#MP = MPRester()
        
configs = configs_reader('config.yaml')

id = "mp-4226"#"Enter the materials project id e.g. mp-4226: "


job = "#!/bin/bash\n" \
      + "#SBATCH --job-name=r" + id + "\n" \
      + "#SBATCH --account=" + configs[1]['account'] + "\n" \
      + "#SBATCH --partition="+configs[1]['partition']+"\n" \
      + "#SBATCH --time="+str(configs[1]['time'])+"\n" \
      + "#SBATCH --ntasks=" + str(configs[1]['cores']) + "\n" \
      + "#SBATCH --mem-per-cpu="+configs[1]['mem-per-cpu']+"\n" \
      + "\n" \
      + "\n" \
      + "module load vasp\n" \
      + "srun vasp_std"


# Make input folders
current_dir = os.getcwd()
os.mkdir(current_dir+ "/relaxation")
os.chdir(current_dir+ "/relaxation")
# 1.1 get conventional unit cell from MP and make POSCAR-unitcell
MP = MPRester(api_key=configs[0])
getbandgap = MP.materials.electronic_structure.search(material_ids=id,fields=['band_gap','magnetic_ordering'])
magnetic_ordering = getbandgap[0].magnetic_ordering
egap = getbandgap[0].band_gap

getp = MP.get_structure_by_material_id(id,conventional_unit_cell=True)
ps = Poscar(getp)
ps.write_file('POSCAR',significant_figures=16)

struct = Structure.from_file('POSCAR')

MPinfo = MPRelaxSet(structure=struct)
MPinfo.potcar.write_file('POTCAR')
MP_INCAR_sett = MPinfo.incar.as_dict()
MAGMOM = MP_INCAR_sett['MAGMOM']

magmom = ''
for i in MAGMOM:
    magmom = magmom + ' ' +str(i)


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


myincar = "PREC = Accurate\n" \
              + "GGA = PS\n" \
              + "IBRION = 2\n" \
              + "NSW = 100\n" \
              + "NELMIN = 5\n" \
              + "ISIF = 3\n" \
              + "ENCUT = " + str(format(1.3 * ENMAX, ".6f")) + "\n" \
              + "EDIFF = 1.000000e-08\n" \
              + "EDIFFG = -1.000000e-06\n" \
              + "IALGO = 38\n" \
              + "LREAL = .FALSE.\n" \
              + "ADDGRID = .TRUE.\n" \
              + "LWAVE = .FALSE.\n" \
              + "LCHARG = .FALSE.\n" \
              + "NPAR = " + str(int(configs[1]['cores'] / 10)) + "\n"


if magnetic_ordering == 'NM' and egap !=0:
    print("Egap ="+ str(egap) +"eV")
    myincar = myincar + "ISMEAR = 0\n" + "SIGMA = 1.000000e-02\n"

elif magnetic_ordering == 'NM' and egap ==0:
    print("Material is metalic")
    myincar = myincar + "ISMEAR = 1\n" + "SIGMA = 0.2\n"

elif magnetic_ordering != "NM" and egap !=0:
    print("The magnetic ordering of this materiasl is: "+ magnetic_ordering)
    print("Egap =" + str(egap) + "eV")
    myincar = myincar + "ISMEAR = 0\n" \
              + "SIGMA = 1.000000e-02\n" \
              + "ISPIN = 2.\n" \
              + "MAGMOM = " + magmom+ "\n"

elif magnetic_ordering != "NM" and egap ==0:
    print("Material is metalic")
    print("The magnetic ordering of this materiasl is: "+ magnetic_ordering)
    myincar = myincar + "ISMEAR = 1\n" \
              + "SIGMA = 0.2\n" \
              + "ISPIN = 2.\n" \
              + "MAGMOM = " + magmom+ "\n"



open('INCAR', 'wt').write(myincar)


# 2.2 KPOINTS
if egap > 2:
    Rk=20
elif egap < 1:
    Rk=40
else:
    Rk=30

N1 = int(max(1,Rk*struct.lattice.reciprocal_lattice.a/(2*np.pi)+0.5))
N2 = int(max(1,Rk*struct.lattice.reciprocal_lattice.b/(2*np.pi)+0.5))
N3 = int(max(1,Rk*struct.lattice.reciprocal_lattice.c/(2*np.pi)+0.5))
mykpt = 'Automatic mesh\n'+'0\n'+'Gamma\n'+ str(N1)+" "+str(N2)+" "+str(N3)+'\n'+'0 0 0'
fout = open("KPOINTS", 'wt').write(mykpt)



open('job_relax.sh', 'wt').write(job)


#Run Calculations
print('Calculation is running...')
subprocess.run('sbatch --wait job_relax.sh', shell=True)

convergence = False
jos_res = 0
while force_check("OUTCAR") == False or convergence == False:
    if not os.path.isfile("POSCAR.org"):
        #print("backup to POSCAR.org")
        os.rename("POSCAR", "POSCAR.org")

    shutil.copy("CONTCAR", "POSCAR")

    for file in ["OUTCAR", "CONTCAR", "XDATCAR", "vasprun.xml"]:
        for i in range(1, 10):
            if os.path.isfile(f"{file}.{i}"):
                print(f"{file}.{i} exists")
            else:
                #print(f"backup to {file}.{i}")
                os.rename(file, f"{file}.{i}")
                break

    jos_res = jos_res + 1
    print("calculation restarting number: "+str(jos_res))
    subprocess.run('sbatch --wait job_relax.sh', shell=True)
    # check Ionic and electronic convergence
    try:
       vrun = Vasprun(filename="vasprun.xml")
       if vrun.converged == True:
            convergence = True
            
    except:
        convergence = False
            





print("Optimization Finished!")
