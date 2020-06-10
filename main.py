import pathlib
import os
import json
from datetime import datetime
from Methods.create_config_file import config_gen
from Methods.utils import launch_sdp, open_pdb_file, env_set
from Methods.distance_file_gen import gen_distance_file
from Methods.spg import launch_spg
from Methods.statistics import launch_statistics
import glob
import logging


def main():
    """This code aims the initial start configuration of the code as also creates an statistical file."""

    # Get the current working dir. path;
    current_dir = str(pathlib.Path().absolute())

    # Get current time
    now = datetime.now()
    current_time = now.strftime("%d-%m-%Y") + '_' + now.strftime("%H_%M_%S")

    # ---------------------------- Set LOG environment ----------------------------
    # set up logging to file - see previous section for more details
    logger_path = f'{current_dir}/LOGs/{current_time}.log'
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-20s: %(levelname)-8s %(message)s',
                        filename=logger_path,
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('[%(name)-20s]: [%(levelname)s] %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Now, we can log to the root logger, or any other logger. First the root...
    logging.info(':: The Process started.')

    # Now, define a couple of other loggers which might represent areas in your
    # application:

    distance = logging.getLogger('main.distance-file')
    # ---------------------------- End LOG environment ----------------------------

    # reads the config file, if it exists;
    config_file_path = current_dir + "\\general_config.txt"
    # if config file exists continue, else create a new one;
    if os.path.isfile(config_file_path):
        # Get the available configuration file;
        with open(config_file_path) as file:
            config = json.load(file)
        logging.info(":: Read configuration file complete.\n")
        pass
    else:
        logging.info(":: File not found! Creating a new configuration...")
        logging.info(
            ":: Please check the Readme file for information regarding the settings configuration standards.\n")
        config = config_gen(current_dir, config_file_path)

    # Get the configuration settings;
    (proteins_path, is_convex_relax, multi_start_phase, mdjeep_source, protein_black_list,
     global_debug_value) = config.values()

    # get the names and directories for all the available proteins; dict(protein_name: its_dir);
    proteins = dict(zip(os.listdir(config["proteins_path"]), glob.glob(f"{proteins_path}\\*")))  # dict

    # ---------------------------- Now starts the tests ----------------------------
    # Create the Node test directory - Node - run date in current directory
    logging.info(':: Creating the test directory...')
    # First we will save the proteins names and paths to be read by the SDP Matlab script;
    # And also, create the test paths for each one, saving it on a dictionary
    protein_tests = {}

    proteins_filepath_dir = current_dir + f"\\Matlab\\proteins.txt"

    with open(proteins_filepath_dir, 'w+') as file:
        for protein, path in proteins.items():
            # We will not load the proteins on the blacklist;
            if protein in protein_black_list:
                continue
            else:
                directory = f'\\Tests\\{current_time}_{protein}'
                test_path = current_dir + directory
                # Create the directory
                os.mkdir(test_path)
                # Saves the directory
                protein_tests[f"{protein}"] = test_path
                logging.info(f":: Directory '{directory}' created.")

                # Writes this information to Matlab usage
                # Scheme: Node -- node_path -- test_path
                file.write(f"{protein},{path},{test_path}\n")

    # Create the distance files for every available proteins
    for node in protein_tests.keys():
        # distance file generator -- generates the distance file, if it doesn't exists, and returns it's directory:
        try:
            # TODO: Insert the overwrite=distance_overwrite option
            gen_distance_file(node, proteins[f"{node}"])
            distance.info(":: Process completed successfully, waiting for data to be read...\n")

        except OSError as err:
            distance.warning(f":: Distance file generator found an error with node: {node} \n"
                             f":: {err}.")
            distance.warning(":: The process was interrupted!")
            continue

    # ------------ Matlab (YALMIP): SDP Program

    # SDP launch and start phase:
    logging.info(':: Start [SDP] program phase.')
    # prepare a structure file for he matlab bash script (SDP Execution);
    launch_sdp(current_dir)
    logging.info("SPG environment configuration set.\n")

    # ------------ Python Refinement: SPG Program
    num_nodes = len(proteins.keys())
    nd_counter = 1  # node counter
    for node in protein_tests.keys():
        logging.info(f":: #{nd_counter} Node: {node} of {num_nodes}")
        # --- Open PDB File
        node_path = proteins[f'{node}']
        pdb_path = node_path + f'\\{node}.txt'
        dist_path = node_path + "\\dist.txt"
        test_path = protein_tests[f'{node}']
        pdb = open_pdb_file(pdb_path)
        logging.debug(':: Environmental properties successful loaded.')
        distancias, u, v, lb, ub, prop_dist = env_set(dist_path)

        # Now begins the refinement process of each protein
        # SPG launch and start phase
        logging.info(":: Start [SPG] program phase.")
        ops = (prop_dist, is_convex_relax, multi_start_phase, global_debug_value)
        launch_spg(node, pdb, test_path, distancias, lb, ub, u, v, ops)
        # -----
        nd_counter += 1

    # ------------ General statistics: static
    launch_statistics(f'{current_dir}/Statistics/{current_time}.txt', protein_tests)
    # End


if __name__ == "__main__":
    main()
