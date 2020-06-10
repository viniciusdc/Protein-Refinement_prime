import argparse
from Methods.distance_file_gen import gen_distance_file
from Methods.atoms_reordination import atoms_re_ordination
from Methods.spectral_projected_gradient import protein_spg
from Methods.pdf_file_gen import write_pdb_file
from Methods.utils import *
from Methods.obj import *
from LOG import os_display_call
from time import time

parser = argparse.ArgumentParser(description="A protein refinement method...")
parser.add_argument("filename", type=str, help="Input protein name")
parser.add_argument(
    "--distance_overwrite",
    default=False,
    help="Enable the overwrite option for the distance file generator function (default: False)",
)
parser.add_argument(
    "--multi_start",
    default=False,
    help="Enable the multi start option. (default: False)",
)
parser.add_argument(
    "--convex_relax",
    default=False,
    help="Indicates whether the starting point originates from convex relaxation (default: False)",
)
parser.add_argument(
    "--local_dir",
    default="C:\\Users\\viniv\\Documents\\protein_tests\\Nodes",
    help="main directory from Testes/Nodes folder (default: .\\Documents\\protein_tests\\Nodes)",
)
parser.add_argument(
    "--debug_mode",
    default=False,
    help="Enables the option for debug mode in the spectral gradient method",
)
args = parser.parse_args()

# protein filename:
filename = args.filename
# local directory:
local_dir = args.local_dir
# main PDB file directory:
dir_pdb = local_dir + f"\\Teste {filename}\\{filename}.txt"
# Start LOG list:
log = []

try:
    # open PDB file (as an np-array):
    pdb = np.genfromtxt(dir_pdb, dtype="str")
    print(":: PDB file read complete! Initiating distance file generation")
    # num_atom_init :: Initial number of atoms
    num_atom_init = int(len(pdb[:, 1]))

    # log archive
    log.append(":: PDB file read complete! Initiating distance file generation")

except FileNotFoundError:
    print(":: PDB file not found!")
    pdb, num_atom_init = [], 0
    print(":: The process was interrupted")

    # log archive
    log.append(":: PDB file not found!")
    log.append(":: The process was interrupted")
    exit()

# distance file generator -- generates the distance file, if it doesn't exists, and returns it's directory:
try:
    raid_d = gen_distance_file(
        pdb, filename, local_dir, overwrite=args.distance_overwrite
    )

    print(":: Process completed successfully, waiting for data to be read...")
    # open distance file (as an np-array):
    distancias = np.genfromtxt(raid_d, dtype="str")
    print(f":: distance file dist_{filename} read complete!")

    # index vectors (for atom pairs [u, v]):
    u, v = (
        np.array(distancias[:, 0], dtype="int"),
        np.array(distancias[:, 1], dtype="int"),
    )
    # it starts from zero on python
    u = u - np.ones(len(u), dtype="int")
    v = v - np.ones(len(v), dtype="int")

    # lower and upper bounds vectors:
    lb, ub = (
        np.array(distancias[:, 8], dtype="float"),
        np.array(distancias[:, 9], dtype="float"),
    )

    prop_dist = 0
    for k in range(len(u)):
        if int(distancias[k][-1]) != 0:
            prop_dist += 1
    prop_dist = int(prop_dist)

    # Log archive
    log.append(":: Process completed successfully, waiting for data to be read...")
    log.append(f":: distance file dist_{filename} read complete!")

except OSError as err:
    print(f":: Distance file generator found an error {err}")
    print(":: The process was interrupted")
    distancias, u, v, lb, ub, prop_dist = [], [], [], [], [], 0

    # log archive
    log.append(f":: Distance file generator found an error {err}")
    log.append(":: The process was interrupted")
    exit()

# Adjust variables:
Noise, TOL, N, M, w = 1e-1, 1e-6, 2000, 15, np.ones(len(lb))
# Noise :: Degree of disturbance of the expected solution;
# TOL :: tolerance for the SPG;
# N :: maximum accepted number of iterations;
# M :: non monotone parameter of GLL line search;
# w :: weight vector

# Pre-initialization:
# Choice/re-ordination of the available atoms in accordance with the distance data file
solution, total_atoms_ord, string_ord, atoms = atoms_re_ordination(pdb, distancias)
solution = centralizar(solution)
log.append(string_ord)

# Initial point origin:
if args.convex_relax:
    print(":: Initial point --Originated from convex relaxation")
    log.append(":: Initial point --Originated from convex relaxation")
    try:
        # directory to initial point file (convex relax solution)
        raid = local_dir + f"\\Teste {filename}\\relax_scaled_{filename}.txt"
        xi = np.genfromtxt(raid)
        xi = centralizar(xi)
        print(f":: relax_scaled_{filename} file read complete.")
        log.append(f":: relax_scaled_{filename} file read complete.")

        # directory to initial point file (non scaled relax solution)
        raid = local_dir + f"\\Teste {filename}\\relax_non_scaled_{filename}.txt"
        non_scaled = np.genfromtxt(raid)
        non_scaled = centralizar(non_scaled)
        print(f":: relax_non_scaled_{filename} file read complete.")
        log.append(f":: relax_non_scaled_{filename} file read complete.")

        # Starting point from convex relaxation:
        # calculating the objective function value for the non scaled initial point (stress)
        dist_non_scaled = dist_matrix_projection(len(u), u, v, lb, ub, non_scaled)
        fo_non_scaled = stress(non_scaled, dist_non_scaled, u, v, w)

        # calculating the objective function value for the scaled initial point (stress)
        dist_non_scaled = dist_matrix_projection(len(u), u, v, lb, ub, xi)
        fo_scaled = stress(xi, dist_non_scaled, u, v, w)

    except FileNotFoundError:
        print(
            f":: File relax_scaled_{filename} or relax_non_scaled_{filename} not found --convex_relax: False"
        )
        log.append(
            f":: File relax_scaled_{filename} or relax_non_scaled_{filename} not found --convex_relax: False"
        )
        args.convex_relax = False
        fo_scaled, fo_non_scaled, xi = 0, 0, []
else:
    print(":: Generating initial point file using a perturbed expected solution.")
    log.append(":: Generating initial point file using a perturbed expected solution.")
    xi = Noise * np.asarray([np.random.normal(0, 1, len(solution)) for i in range(3)]).T
    xi = xi + solution
    xi = np.array(xi)
    xi = centralizar(xi)
    fo_non_scaled, fo_scaled = 0, 0

# Initiate solver:
if args.multi_start:
    data = {}
    try:
        print(f":: Multi-start option --Enable {10}x times")
        log.append(f":: Multi-start option --Enable {10}x times")
        print(f":: maximum iterations: {N}, tol: {TOL} and memory: {M}")
        log.append(f":: maximum iterations: {N}, tol: {TOL} and memory: {M}")
        # multi start option enable, initial distance vector will be set as in multi-start definition
        yi = dist_matrix_projection(int(len(u)), u, v, lb, ub, xi)
        to = time()
        out = protein_spg(
            stress,
            grad_stress,
            xi,
            yi,
            [u, v, w, lb, ub, TOL, N, M],
            debug_mode=args.debug_mode,
        )
        elapsed_time = time() - to
        fo = stress(xi, yi, u, v, w)

        if args.debug_mode:
            X = out[0]
            if check_solution_dimension(X, solution):
                print(
                    ":: Solution found and expected solution have different number of atoms!"
                )
                log.append(
                    ":: Solution found and expected solution have different number of atoms!"
                )

        print(f":: ({0}) Solution found!")
        data[0] = (out, elapsed_time, fo)

        for i in range(1, 10):
            yi = dist_matrix_projection(int(len(u)), u, v, lb, ub, xi, multistart=True)
            to = time()
            out = protein_spg(
                stress,
                grad_stress,
                xi,
                yi,
                [u, v, w, lb, ub, TOL, N, M],
                debug_mode=args.debug_mode,
            )
            elapsed_time = time() - to
            fo = stress(xi, yi, u, v, w)

            if args.debug_mode:
                X = out[0]
                if check_solution_dimension(X, solution):
                    print(
                        ":: Solution found and expected solution have different number of atoms!"
                    )
                    log.append(
                        ":: Solution found and expected solution have different number of atoms!"
                    )
            print(f":: ({i}) Solution found!")
            data[i] = (out, elapsed_time, fo)

    except Exception as ee:
        print(
            f":: Attempt for protein_spg --Multi-start --Enable failed with --bad error: {ee}"
        )
        print(":: The process was interrupted!")
        log.append(
            f":: Attempt for protein_spg --Multi-start --Enable failed with --bad error: {ee}"
        )
        log.append(":: The process was interrupted!")
        exit()
        out = []
        fo, elapsed_time = 0.0, 0.0

    # parameter initialization for data output:
    ops = (xi, solution, u, v, lb, ub)
    main = (
        filename,
        num_atom_init,
        total_atoms_ord,
        len(u),
        prop_dist,
        args.convex_relax,
        fo_non_scaled,
        fo_scaled,
        ops,
    )
    os_display_call(log, local_dir, main, data, multistart=True)

else:
    # multi start option not enable, initial distance vector will be the standard
    yi = dist_matrix_projection(int(len(u)), u, v, lb, ub, xi)
    if args.debug_mode:
        print(":: Debug mode --True")
        log.append(":: Debug mode --True")
    try:
        print(f":: maximum iterations: {N}, tol: {TOL} and memory: {M}")
        log.append(f":: maximum iterations: {N}, tol: {TOL} and memory: {M}")
        to = time()
        out = protein_spg(
            stress,
            grad_stress,
            xi,
            yi,
            [u, v, w, lb, ub, TOL, N, M],
            debug_mode=args.debug_mode,
        )
        elapsed_time = time() - to
        fo = stress(xi, yi, u, v, w)
        log.append(out[-1])
        print(":: Solution found!")
        log.append(":: Solution found!")

        if args.debug_mode:
            X = out[0]
            if check_solution_dimension(X, solution):
                print(
                    ":: Solution found and expected solution have different number of atoms!"
                )
                log.append(
                    ":: Solution found and expected solution have different number of atoms!"
                )
        data = (out, elapsed_time, fo)

    except Exception as err:
        print(f":: Attempt for protein_spg failed with --bad error: {err}")
        print(":: The process was interrupted!")
        log.append(f":: Attempt for protein_spg failed with --bad error: {err}")
        log.append(":: The process was interrupted!")
        exit()
        out = []
        fo, elapsed_time = 0.0, 0.0
        data = (out, elapsed_time, fo)

    ops = (xi, solution, u, v, lb, ub)
    main = (
        filename,
        num_atom_init,
        total_atoms_ord,
        len(u),
        prop_dist,
        args.convex_relax,
        fo_non_scaled,
        fo_scaled,
        ops
    )
    os_display_call(log, local_dir, main, data)

# ----------------------------------------------------------
#                 Output the data to a file
# ----------------------------------------------------------
atoms = np.array(atoms)


def raid_gen(arch):
    return local_dir + f"\\Teste {filename}\\{arch}_{filename}.pdb"


# export solution pdb file:
# output PDB file format -- ATOM i atom amino A res x1 x2 x3 0.00  0.00 atom[0]:
write_pdb_file(raid_gen("Sol"), out[0], atoms[:, 1], atoms[:, 3], atoms[:, 2])

# export original solution pdb file -- only if it dos not already exists
try:
    with open(raid_gen("Orig")) as f:
        print(f"Orig_{filename} already exists!")
except FileNotFoundError:
    print(f"Creating a new Orig_{filename}")
    write_pdb_file(raid_gen("Orig"), solution, atoms[:, 1], atoms[:, 3], atoms[:, 2])

# export initial point
write_pdb_file(raid_gen("Ponto"), xi, atoms[:, 1], atoms[:, 3], atoms[:, 2])
# ----------------------------------------------------------
