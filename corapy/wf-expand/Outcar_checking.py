from pymatgen.io.vasp.outputs import Outcar
import os
import numpy as np
import shutil

def force_check(inputfile):
    #current_dir = os.getcwd()
    outcar = Outcar(inputfile)

    forces = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        footer_pattern=r"\s--+",
        postprocess=lambda x: float(x),
        last_one_only=False
    )

    force_converge_threshold = 1e-4
    # last step forces
    TOTAL_FORCE = forces[len(forces) - 1]
    all_forces = []
    for i in TOTAL_FORCE:
        for value in i:
            all_forces.append(value)

    #print(all_forces)



    if all(np.abs(x) < force_converge_threshold for x in i) and len(outcar.run_stats) > 2:
        # system relaxed
        return True
    else:
        return False


