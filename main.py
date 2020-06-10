import pathlib
import os
import numpy as np
import json
from Methods.create_config_file import config_gen
import glob


def main():
    """This code aims the initial start configuration of the code as also creates an statistical file."""

    # Get the current working dir. path;
    current_dir = str(pathlib.Path().absolute())

    # reads the config file, if it exists;
    config_file_path = current_dir + "\\general_config.txt"
    try:
        # Check for available configuration file;
        with open(config_file_path) as file:
            config = json.load(file)
    except FileNotFoundError as ee:
        # else, creates a new one;
        print(":: File not found! Creating a new configuration...")
        print(":: Please check the Readme file for information regarding the settings configuration standards.\n")
        config = config_gen(current_dir, config_file_path)

    # Get the configuration settings;
    (proteins_path, is_convex_relax, multi_start_phase, mdjeep_source, protein_black_list, global_debug_value) = config

    # get the names and directories for all the available proteins; dict(protein_name: its_dir);
    proteins = dict(zip(os.listdir(config["proteins_path"]), glob.glob(f"{proteins_path}\\*")))  # dict

    # Now starts the tests:
    for protein, path in proteins.items():


# prepare a structure file for he matlab bash script (SDP Execution);


if __name__ == "__main__":
    main()
