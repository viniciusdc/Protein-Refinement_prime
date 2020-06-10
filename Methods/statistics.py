import logging
import json


def launch_statistics(out: str, tests: dict) -> None:
    """This function unite all the available information from the tested nodes and outputs into a file"""
    # gen statistics from testes path

    with open(out, 'w+') as static:
        static.write("#######################################  Statistics  #######################################\n")
        for Node, node_test_path in tests.items():
            # Get PDB results from spg-static-multi-start-LOG file
            spg_static_path = node_test_path + '\\spg_static_multistart_LOG.txt'
            sdp_static_path = node_test_path + '\\solver_varargout.txt'
            with open(spg_static_path) as file:
                spg_static = json.load(file)
            with open(sdp_static_path) as file:
                sdp_static = json.load(file)

            num_atom_init = spg_static["init_atom_#"]
            total_atoms_ord = spg_static["atom_#_re-ordination"]
            m = spg_static["assessed_dist"]
            prop_dist = spg_static["Know_dist"]
            static.write(
                f":: Protein: {Node}, Initial atoms number: {num_atom_init}, after re-ordination {total_atoms_ord}.\n"
            )
            static.write(f":: Assessed distances: {m} and known distances: {prop_dist}.\n")

            convex = spg_static["convex"]
            if convex:
                static.write(":: [SPD] Results:\n")
                fo_non_scaled = spg_static["init_fun_val_relax"]
                fo_scaled = spg_static["init_fun_val_relax_k"]
                static.write(
                    f":: Initial objective value for the relaxed problem: {fo_non_scaled}.\n"
                )
                static.write(
                    f":: Initial objective value for the relaxed problem --scaled {fo_scaled}.\n"
                )
                static.write(f":: Yalmip Time (s): {sdp_static['yalmiptime']:<12} "
                             f"Solvertime (Sedumi): {sdp_static[f'solvertime']:<12} "
                             f"Total time (s): {sdp_static[f'elapsed time']}\n")

            rmsd_i = spg_static["RMSDi"]
            mde_i = spg_static["MDEi"]
            static.write(f":: RMSDi = {rmsd_i:<24} MDEi = {mde_i}\n")

            multi_start = spg_static["multi-start"]
            if type(multi_start) == list:
                static.write(":: spg results --multi start: True\n")
                static.write(
                    ":: Iter - bck -- RMSDf ----- MDEf"
                    " ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)\n"
                )
                for attempts in multi_start:
                    it, back, rmsdf, mdef, fun_i, fun_f, gtd, norm_d, time = attempts.values()
                    static.write(
                        f"   {it}: {back} {rmsdf:<11} {mdef:<10} {fun_i:<11} "
                        f"{fun_f:<10} {gtd:<10} {norm_d:<10} {time}\n"
                    )
            elif type(spg_static["standard"]) == dict:
                static.write(":: spg results --multi start: False\n")

                static.write(
                    ":: Iter - bck -- RMSDf ----- MDEf"
                    " ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)\n"
                )
                it, back, rmsdf, mdef, fun_i, fun_f, gtd, norm_d, time = spg_static["standard"]
                static.write(
                        f"   {it}: {back} {rmsdf:<11} {mdef:<10} {fun_i:<11} "
                        f"{fun_f:<10} {gtd:<10} {norm_d:<10} {time}\n"
                    )
            static.write(
                "############################################################################################\n"
            )
    return
