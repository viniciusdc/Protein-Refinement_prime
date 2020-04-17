from Methods.utils import rmsd, mde
from datetime import datetime
import sys


def os_display_call(log, local_dir, main, data, multistart=False):
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

    print("#####################################  INFO  #####################################")
    print(
        f":: Protein: {filename}, Initial atoms number: {num_atom_init}, after re-ordination {total_atoms_ord}."
    )
    print(f":: Assessed distances: {m} and known distances: {prop_dist}.")
    if convex:
        print(f":: Initial objective value for the relaxed problem: {fo_non_scaled:.4e}")
        print(f':: Initial objective value for the relaxed problem --scaled {fo_scaled:.4e}')
    rmsd_i, mde_i = rmsd(xi, solution), mde(xi, u, v, lb, ub)
    print(f":: RMSDi = {rmsd_i:<24.4e} MDEi = {mde_i:.4e}")

    # -----------------------------------------------------------------------------------
    # Multi-start option --Enabled
    # -----------------------------------------------------------------------------------
    if multistart:
        if type(data) != dict:
            print(":: data type object not match with dict structure!")
            print(":: The process was interrupted")
            return exit()
        print(":: spg results --multi start: True")
        print(":: Iter -- bck -- RMSDf ----- MDEf"
              " ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)")
        sub_log = {}
        for key in data:
            out, elapsed_time, fo = data[key]
            x_spg, backtracking, iterations, fun_o, gtd, norm_d, log_spg = out
            # Statistics:
            rmsd_f = rmsd(x_spg, solution)
            mde_f = mde(x_spg, u, v, lb, ub)
            prompt_string = f'   {iterations:<7} {backtracking:<6} {rmsd_f:<11.2e} {mde_f:<10.2e} {fo / 2:<11.2e} ' \
                            f'{fun_o / 2:<11.2e} {gtd:<10.2e} {norm_d:<10.2e} {elapsed_time:.3f}'
            sub_log[key] = prompt_string
            print(prompt_string)
        print("##################################################################################")

        # -----------------------------------------------------------------------------
        # Generating output file with statistics:
        # -----------------------------------------------------------------------------
        now = datetime.now()
        dt_time = now.strftime('%d-%m-%Y--%H-%M-%S')
        out_file_log = local_dir + f'\\Teste {filename}\\LOG_({dt_time}).txt'
        with open(out_file_log, 'w') as f:
            for item in log:
                if type(item) == list:
                    for item_spg in item:
                        f.write(f"{item_spg}\n")
                    continue
                f.write(f"{item}\n")
            f.write("#####################################  INFO  #####################################\n")
            f.write(
                f":: Protein: {filename}, Initial atoms number: {num_atom_init}, "
                f"after re-ordination {total_atoms_ord}.\n"
            )
            f.write(
                f":: Assessed distances: {m} and known distances: {prop_dist}.\n"
            )
            if convex:
                f.write(f":: Initial objective value for the relaxed problem: {fo_non_scaled:.4e}\n")
                f.write(f':: Initial objective value for the relaxed problem --scaled {fo_scaled:.4e}\n')
            f.write(f":: RMSDi = {rmsd_i:<24.4e} MDEi = {mde_i:.4e}\n")
            if type(data) != dict:
                print(":: data type object not match with dict structure!\n")
                print(":: The process was interrupted\n")
            f.write(":: spg results --multi start: True\n")
            f.write(":: Iter -- bck -- RMSDf ----- MDEf ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)\n")
            for key in sub_log:
                f.write(f'{sub_log[key]}\n')
            f.write("##################################################################################")

    # -----------------------------------------------------------------------------------
    # Multi-start --Disable  Standard
    # -----------------------------------------------------------------------------------
    else:
        out, elapsed_time, fo = data
        x_spg, backtracking, iterations, fun_o, gtd, norm_d, log_spg = out
        # Statistics:
        rmsd_f = rmsd(x_spg, solution)
        mde_f = mde(x_spg, u, v, lb, ub)
        print(":: spg results --multi start: False")
        print(":: Iter -- bck -- RMSDf ----- MDEf ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)")
        prompt_string = f'   {iterations:<7} {backtracking:<6} {rmsd_f:<11.2e} {mde_f:<10.2e} {fo / 2:<11.2e} ' \
                        f'{fun_o / 2:<11.2e} {gtd:<10.2e} {norm_d:<10.2e} {elapsed_time:.3f}'
        print(prompt_string)
        print("##################################################################################")

        # -----------------------------------------------------------------------------
        # Generating output file with statistics:
        # -----------------------------------------------------------------------------
        now = datetime.now()
        dt_time = now.strftime('%d-%m-%Y--%H-%M-%S')
        out_file_log = local_dir + f'\\Teste {filename}\\LOG_({dt_time}).txt'
        with open(out_file_log, 'w') as f:
            for item in log:
                if type(item) == list:
                    for item_spg in item:
                        f.write(f"{item_spg}\n")
                    continue
                f.write(f"{item}\n")

            f.write("#####################################  INFO  #####################################\n")
            f.write(
                f":: Protein: {filename}, Initial atoms number: {num_atom_init}, "
                f"after re-ordination {total_atoms_ord}.\n"
            )
            f.write(
                f":: Assessed distances: {m} and known distances: {prop_dist}.\n"
            )
            if convex:
                f.write(f":: Initial objective value for the relaxed problem: {fo_non_scaled:.3e}\n")
                f.write(f':: Initial objective value for the relaxed problem --scaled {fo_scaled:.3e}\n')
            else:
                f.write(f":: Initial objective value: {fo:.3e}\n")
            f.write(f":: RMSDi = {rmsd_i:<24.3e} MDEi = {mde_i:.3e}\n")
            f.write(":: spg results --multi start: False\n")
            # statistics
            f.write(":: Iter -- bck -- RMSDf ----- MDEf ----- i_val ----- f_val ----- gtd ----- |d| ----- time(s)\n")
            f.write(f'{prompt_string}\n')
            f.write("##################################################################################\n")
