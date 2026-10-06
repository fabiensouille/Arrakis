#!/usr/bin/env python3
"""
Run 

How to execute: arrakis.py <my_case>.yml
"""
from pathlib import Path
import os
import platform
import subprocess
import argparse
import time
import glob
import yaml

ARRAKIS_DIR = Path(os.path.abspath(__file__)).parents[2]
SCR_DIR = os.path.join(ARRAKIS_DIR, "scripts", "arrakis")
SRC_DIR = os.path.join(ARRAKIS_DIR, "sources")
BIN_DIR = os.path.join(ARRAKIS_DIR, "build", "bin")

# IO functions for cas.yml input file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_yml_input_file(file):
    """
    Read case parameters file
    """
    with open(file, "r") as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            return 0

def write_yml_input_file(file, case_dict, version='2.0'):
    """
    Write case parameters file
    """
    f = open(file, "a")
    f.write("\n")
    f.write("\n")
    f.write("#====================================================#\n")
    f.write("#  CASE INPUT DATA <ARRAKIS> ${}                     #\n".format(version))
    f.write("#====================================================#\n")
    f.write("\n")

    for section, secdict in case_dict.items():
        f.write("\n")
        f.write("#====================================================#\n")
        f.write("{}:\n".format(section))
        f.write("#====================================================#\n")
        for key, value in secdict.items():
            f.write("  {:<51}: {}\n".format(key, value))
        f.write("\n")
    f.close()

def update_dico(default_dict, custom_dict):
    """
    Update dictionnary without removed existing keys.
    """
    # update custom_dict with default values
    for item in default_dict:
        for value in default_dict[item].keys():
            if value not in custom_dict[item].keys():
                custom_dict[item].update({value: default_dict[item][value]})
    # sort custom dict
    output_dict = default_dict
    for item in output_dict:
        for value in output_dict[item].keys():
            output_dict[item].update({value: custom_dict[item][value]})
    return output_dict

# Run functions
# ~~~~~~~~~~~~~
def run_arrakis(case_dict, omp_num_threads=4, profiling=True):
    """
    Run ARRAKIS
    """
    # read refault yaml input file
    case_inputs = read_yml_input_file(os.path.join(SCR_DIR, "arrakis.yml"))
    case_inputs = update_dico(case_inputs, case_dict)

    # remove existing run.yml if present, then write a new one
    Path("run.yml").unlink(missing_ok=True)
    write_yml_input_file("run.yml", case_inputs)

    # create directories
    os.makedirs("RESU", exist_ok=True)
    os.makedirs("FIG", exist_ok=True)

    # run ARRAKIS
    time_ini = time.time()
    # Handle Windows executable extension
    arrakis_exe = "ARRAKIS.exe" if platform.system() == "Windows" else "ARRAKIS"
    arrakis_path = os.path.join(BIN_DIR, arrakis_exe)
    
    # Set environment variables and run
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(omp_num_threads)
    subprocess.run([arrakis_path], env=env, check=True)
    time_fin = time.time()

    if profiling:
        print("CPU Time = ", time_fin - time_ini)

    # remove run.yml
    Path("run.yml").unlink(missing_ok=True)

if __name__ == '__main__':
    # parse options
    parser = argparse.ArgumentParser(description='Run ARRAKIS')
    parser.add_argument('case', type=str, nargs='+', help='input arrakis data file (.yml)')
    parser.add_argument('-r', '--run', default=False,
                        action="store_true", help='run arrakis')
    parser.add_argument('-c', '--clean', default=False,
                        action="store_true", help='clean RESU dir')
    parser.add_argument('-nc', '--ncsize', type=int, default=4,
                        help='number of OpenMP threads (default: 4)')

    ARGS = parser.parse_args()
    NO_ARGS = bool(ARGS.run is False and ARGS.clean is False)

    # read yaml input file
    CASE_DICT = read_yml_input_file(ARGS.case[0])

    # run ARRAKIS
    if NO_ARGS or ARGS.run:
        run_arrakis(CASE_DICT, omp_num_threads=ARGS.ncsize)

    # post processing
    if ARGS.clean:
        for file in glob.glob("RESU/*"):
            if os.path.isfile(file):
                os.remove(file)
