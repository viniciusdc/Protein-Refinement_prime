import logging
import numpy as np
from Methods.atoms_reordination import atoms_re_ordination
from Methods.spectral_projected_gradient import protein_spg
from Methods.pdf_file_gen import write_pdb_file
from Methods.utils import *
from Methods.obj import *
from LOG import os_display_call
from time import time


def launch_spg(node, pdb, test_path, distancias, lb, ub, u, v, config_ops):
    """This method aims to create all the environment variables to execute the spectral gradient method."""
    # Get run options
    (prop_dist, is_convex_relax, multi_start, debug_mode) = config_ops

    # Get logger
    logger = logging.getLogger('root.load_spg')
    if debug_mode:
        logger.setLevel(10)

    # Adjust variables:
    Noise, TOL, N, M, w = 1e-1, 1e-6, 2000, 15, np.ones(len(lb))
    num_atom_init = int(len(pdb[:, 1]))
    # Noise :: Degree of disturbance of the expected solution;
    # TOL :: tolerance for the SPG;
    # N :: maximum accepted number of iterations;
    # M :: non monotone parameter of GLL line search;
    # w :: weight vector

    # Pre-initialization:
    # Choice/re-ordination of the available atoms in accordance with the distance data file
    solution, total_atoms_ord, atoms = atoms_re_ordination(pdb, distancias)
    solution = centralizar(solution)

    # Initial point origin:
    if is_convex_relax:
        logger.info(":: Initial point --Originated from convex relaxation")
        try:
            # directory to initial point file (convex relax solution)
            raid = test_path + f"\\relax_k.txt"
            xi = np.genfromtxt(raid)
            xi = centralizar(xi)
            logger.debug(f":: relax_k file read complete.")

            # directory to initial point file (non scaled relax solution)
            raid = test_path + f"\\relax_np.txt"
            non_scaled = np.genfromtxt(raid)
            non_scaled = centralizar(non_scaled)
            logger.debug(f":: relax_np file read complete.")

            # Starting point from convex relaxation:
            # calculating the objective function value for the non scaled initial point (stress)
            dist_non_scaled = dist_matrix_projection(len(u), u, v, lb, ub, non_scaled)
            fo_non_scaled = stress(non_scaled, dist_non_scaled, u, v, w)

            # calculating the objective function value for the scaled initial point (stress)
            dist_non_scaled = dist_matrix_projection(len(u), u, v, lb, ub, xi)
            fo_scaled = stress(xi, dist_non_scaled, u, v, w)

        except FileNotFoundError:
            logger.warning(f":: File relax_k or relax_np not found --convex_relax: False")
            is_convex_relax = False
            fo_scaled, fo_non_scaled, xi = 0, 0, []
    else:
        logger.debug(":: Generating initial point file using a perturbed expected solution.")
        xi = Noise * np.asarray([np.random.normal(0, 1, len(solution)) for i in range(3)]).T
        xi = xi + solution
        xi = np.array(xi)
        xi = centralizar(xi)
        fo_non_scaled, fo_scaled = 0, 0

    # Initiate solver:
    if multi_start:
        data = {}
        try:
            logger.info(f":: Multi-start option --Enable {10}x times")
            logger.info(f":: maximum iterations: {N}, tol: {TOL} and memory: {M}")
            # multi start option enable, initial distance vector will be set as in multi-start definition
            yi = dist_matrix_projection(int(len(u)), u, v, lb, ub, xi)
            to = time()
            out = protein_spg(
                stress,
                grad_stress,
                xi,
                yi,
                [u, v, w, lb, ub, TOL, N, M],
                debug_mode=debug_mode,
            )
            elapsed_time = time() - to
            fo = stress(xi, yi, u, v, w)
            X = out[0]
            if check_solution_dimension(X, solution):
                logger.debug(":: Solution found and expected solution have different number of atoms!")

            logger.info(f":: ({0}) Solution found!")
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
                    debug_mode=debug_mode,
                )
                elapsed_time = time() - to
                fo = stress(xi, yi, u, v, w)

                X = out[0]
                if check_solution_dimension(X, solution):
                    logger.debug(
                        ":: Solution found and expected solution have different number of atoms!"
                        )
                logger.info(f":: ({i}) Solution found!")
                data[i] = (out, elapsed_time, fo)
            logger.info(f':: Multi-start phase of node: {node} complete.\n')

        except Exception as ee:
            logger.info(":: The process was interrupted!")
            logger.warning(
                f":: Attempt for protein_spg --Multi-start --Enable failed with --bad error: {ee}"
            )
            exit()
            out = []
            fo, elapsed_time = 0.0, 0.0

        # parameter initialization for data output:
        ops = (xi, solution, u, v, lb, ub)
        main = (
            node,
            num_atom_init,
            total_atoms_ord,
            len(u),
            prop_dist,
            is_convex_relax,
            fo_non_scaled,
            fo_scaled,
            ops,
        )
        os_display_call(test_path, main, data, multistart=True)

    else:
        # multi start option not enable, initial distance vector will be the standard
        yi = dist_matrix_projection(int(len(u)), u, v, lb, ub, xi)
        if debug_mode:
            logger.debug(":: --True")
        try:
            logger.info(f":: maximum iterations: {N}, tol: {TOL} and memory: {M}")
            to = time()
            out = protein_spg(
                stress,
                grad_stress,
                xi,
                yi,
                [u, v, w, lb, ub, TOL, N, M],
                debug_mode=debug_mode,
            )
            elapsed_time = time() - to
            fo = stress(xi, yi, u, v, w)
            logger.info(":: Solution found!")

            if debug_mode:
                X = out[0]
                if check_solution_dimension(X, solution):
                    logger.debug(
                        ":: Solution found and expected solution have different number of atoms!"
                    )
            data = (out, elapsed_time, fo)
            logger.info(f':: Standard phase of node: {node} complete.\n')

        except Exception as err:
            logger.warning(f":: Attempt for protein_spg failed with --bad error: {err}")
            logger.warning(":: The process was interrupted!")
            exit()
            out = []
            fo, elapsed_time = 0.0, 0.0
            data = (out, elapsed_time, fo)

        ops = (xi, solution, u, v, lb, ub)
        main = (
            node,
            num_atom_init,
            total_atoms_ord,
            len(u),
            prop_dist,
            is_convex_relax,
            fo_non_scaled,
            fo_scaled,
            ops
        )
        os_display_call(test_path, main, data)

    # ----------------------------------------------------------
    #                 Output the data to a file
    # ----------------------------------------------------------
    atoms = np.array(atoms)

    def raid_gen(arch):
        return test_path + f"\\{arch}.pdb"

    # export solution pdb file:
    # output PDB file format -- ATOM i atom amino A res x1 x2 x3 0.00  0.00 atom[0]:
    write_pdb_file(raid_gen("Sol"), out[0], atoms[:, 1], atoms[:, 3], atoms[:, 2])

    # export original solution pdb file -- only if it dos not already exists
    try:
        with open(raid_gen("Orig")) as f:
            logger.info(f":: Origin already exists!")
    except FileNotFoundError:
        logger.info(f":: Creating a new Origin\n")
        write_pdb_file(raid_gen("Orig"), solution, atoms[:, 1], atoms[:, 3], atoms[:, 2])

    # export initial point
    write_pdb_file(raid_gen("Ponto"), xi, atoms[:, 1], atoms[:, 3], atoms[:, 2])
    # ----------------------------------------------------------

    return
