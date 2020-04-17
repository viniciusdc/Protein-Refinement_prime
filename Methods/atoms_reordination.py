import numpy as np
from time import time


def atoms_re_ordination(pdb_file, dist_file):
    """sorting atoms according to its distance file, to create the expected/correct solution"""
    # total number of atoms (in accordance with the distance data file)
    total_atoms_ord = max(
        max(np.array(dist_file[:, 0], dtype="int")),
        max(np.array(dist_file[:, 1], dtype="int")),
    )

    to = time()
    atoms_names_list = []
    for item in dist_file:
        line = [item[0]] + list(item[2:5])
        atoms_names_list.append(line)
        line = [item[1]] + list(item[5:8])
        atoms_names_list.append(line)

    atoms = []
    for k in range(1, total_atoms_ord + 1):
        for item in atoms_names_list:
            if k == int(item[0]):
                atoms.append(item)
                break

    # creating correct solution point
    atoms_solve = []
    for item in atoms:
        for atom in pdb_file:
            atom_name = [atom[2], atom[3], atom[5]]
            if all(item[k + 1] == atom_name[k] for k in range(3)):
                atoms_solve.append(np.array(atom[6:9], dtype="float"))

    solution = np.array(atoms_solve)
    elapsed_time = time() - to
    string = f":: The process of reformatting/re-ordination was successfully completed in {elapsed_time:.4f} s"
    print(string)

    return solution, total_atoms_ord, string, atoms
