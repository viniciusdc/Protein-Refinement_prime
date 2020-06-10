from Methods.utils import rmsd, mde
from datetime import datetime
import logging
import json
import sys


def os_display_call(test_path, main, data, multistart=False):
    (
        filename,
        num_atom_init,
        total_atoms_ord,
        m,
        prop_dist,
        convex,
        fo_non_scaled,
        fo_scaled,
        ops,
    ) = main
    xi, solution, u, v, lb, ub = ops
    # Get logger
    logger = logging.getLogger('root.spgLOG')
    logger.info(
        "##########################################  INFO  ##########################################"
    )
    logger.info(
        f":: Protein: {filename}, Initial atoms number: {num_atom_init}, after re-ordination {total_atoms_ord}."
    )
    logger.info(f":: Assessed distances: {m} and known distances: {prop_dist}.")
    if convex:
        logger.info(
            f":: Initial objective value for the relaxed problem: {fo_non_scaled:.4e}"
        )
        logger.info(
            f":: Initial objective value for the relaxed problem --scaled {fo_scaled:.4e}"
        )
    rmsd_i, mde_i = rmsd(xi, solution), mde(xi, u, v, lb, ub)
    logger.info(f":: RMSDi = {rmsd_i:<24.2e} MDEi = {mde_i:.2e}")

    # -----------------------------------------------------------------------------------
    # Multi-start option --Enabled
    # -----------------------------------------------------------------------------------
    if multistart:
        if type(data) != dict:
            logger.warning(":: data type object not match with dict structure!")
            logger.warning(":: The process was interrupted")
            return exit()
        logger.info(":: spg results --multi start: True")
        logger.info(
            ":: Iter - bck -- RMSDf ----- MDEf"
            " ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)"
        )
        sub_log = {}
        k = 0
        for key in data:
            out, elapsed_time, fo = data[key]
            x_spg, backtracking, iterations, fun_o, gtd, norm_d = out
            # Statistics:
            rmsd_f = rmsd(x_spg, solution)
            mde_f = mde(x_spg, u, v, lb, ub)
            prompt_string = (
                f"   {iterations:<5}: {backtracking:<6} {rmsd_f:<11.2e} {mde_f:<10.2e} {fo / 2:<11.2e} "
                f"{fun_o / 2:<10.2e} {gtd:<10.2e} {norm_d:<10.2e} {elapsed_time:.3f}"
            )
            sub_log[k] = {"iter": f'{iterations:<7}', "back": f'{backtracking:<6}', "RMDSf": f'{rmsd_f:<11.2e}',
                          "MDEf": f'{mde_f:<10.2e}', "fun_i": f'{fo / 2:<11.2e}', "fun_f": f'{fun_o / 2:<10.2e}',
                          "gtd": f'{gtd:<10.2e}', "norm_d": f'{norm_d:<10.2e}', "time": f'{elapsed_time:.3f}'}
            logger.info(prompt_string)
            k += 1

        logger.info(
            "############################################################################################"
        )

        # -----------------------------------------------------------------------------
        # Generating output file with statistics:
        # -----------------------------------------------------------------------------
        static_dict = {"node": f'{filename}', "init_atom_#": f"{num_atom_init}",
                       "atom_#_re-ordination": f'{total_atoms_ord}',
                       "assessed_dist": f'{m}', "Know_dist": f'{prop_dist}'}
        if convex:
            static_dict["convex"] = True
            static_dict["init_fun_val_relax"] = f'{fo_non_scaled:.4e}'
            static_dict["init_fun_val_relax_k"] = f'{fo_scaled:.4e}'
        else:
            static_dict["convex"] = False
            static_dict["init_fun_val_relax"] = 'N/A'
            static_dict["init_fun_val_relax_k"] = 'N/A'

        static_dict["RMSDi"] = f'{rmsd_i:<24.2e}'
        static_dict["MDEi"] = f'{mde_i:.2e}'
        if type(data) != dict:
            logger.warning(":: data type object not match with dict structure!\n")
            logger.warning(":: The process was interrupted\n")
        multistart_list = []
        n = len(sub_log.keys())
        for i in range(n):
            multistart_list.append(sub_log[i])
        static_dict["multi-start"] = multistart_list
        static_dict["standard"] = False

        static_log = test_path + f"\\spg_static_multistart_LOG.txt"
        with open(static_log, "w") as f:
            json.dump(static_dict, f)

    # -----------------------------------------------------------------------------------
    # Multi-start --Disable  Standard
    # -----------------------------------------------------------------------------------
    else:
        out, elapsed_time, fo = data
        x_spg, backtracking, iterations, fun_o, gtd, norm_d = out
        # Statistics:
        rmsd_f = rmsd(x_spg, solution)
        mde_f = mde(x_spg, u, v, lb, ub)
        logger.info(":: spg results --multi start: False")
        logger.info(
            ":: Iter - bck -- RMSDf ----- MDEf"
            " ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)"
        )
        prompt_string = (
            f"   {iterations:<5}: {backtracking:<6} {rmsd_f:<11.2e} {mde_f:<10.2e} {fo / 2:<11.2e} "
            f"{fun_o / 2:<10.2e} {gtd:<10.2e} {norm_d:<10.2e} {elapsed_time:.3f}"
        )

        logger.info(prompt_string)
        logger.info(
            "############################################################################################"
        )

        # -----------------------------------------------------------------------------
        # Generating output file with statistics:
        # -----------------------------------------------------------------------------
        static_log = test_path + f"\\spg_static_standard_LOG.txt"

        static_dict = {"node": f'{filename}', "init_atom_#": f"{num_atom_init}",
                       "atom_#_re-ordination": f'{total_atoms_ord}',
                       "assessed_dist": f'{m}', "Know_dist": f'{prop_dist}'}
        if convex:
            static_dict["convex"] = True
            static_dict["init_fun_val_relax"] = f'{fo_non_scaled:.4e}'
            static_dict["init_fun_val_relax_k"] = f'{fo_scaled:.4e}'
        else:
            static_dict["convex"] = False
            static_dict["init_fun_val_relax"] = 'N/A'
            static_dict["init_fun_val_relax_k"] = 'N/A'

        static_dict["RMSDi"] = f'{rmsd_i:<24.2e}'
        static_dict["MDEi"] = f'{mde_i:.2e}'
        static_dict["multi-start"] = False
        static_dict["standard"] = {"iter": f'{iterations:<7}', "back": f'{backtracking:<6}',
                                   "RMDSf": f'{rmsd_f:<11.2e}',
                                   "MDEf": f'{mde_f:<10.2e}', "fun_i": f'{fo / 2:<11.2e}',
                                   "fun_f": f'{fun_o / 2:<10.2e}',
                                   "gtd": f'{gtd:<10.2e}', "norm_d": f'{norm_d:<10.2e}',
                                   "time": f'{elapsed_time:.3f}'}

        with open(static_log, "w") as file:
            json.dump(static_dict, file)
