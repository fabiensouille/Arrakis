#!/usr/bin/env python3
"""
Run validation

How to execute: validate.py
"""
import os
import sys
import glob
import subprocess
import argparse
from pathlib import Path

ARRAKIS_DIR = Path(os.path.abspath(__file__)).parents[2]

def validate_sv(ini=False):
    """
    Run test located in examples
    """
    caselist = [\
        "SV_00_dambreak",
        "SV_01_bumpsub",
        "SV_02_bumpcri",
        "SV_03_bumptrans",
        "SV_04_thacker",
        "SV_05_lakeatrest",
        "SV_07_flume_without_friction",
        "SV_06_gouttedeau",
        "SV_08_flume_with_friction",
        "SV_09_flume_periodic",
        "SV_10_tracer",
        "SV_11_tracer_multiple",
        "SV_12_travelling_jump",
        "SV_13_perturbations",
        "SV_14_solitary_wave",
        "SV_15_run_up",
        "SV_16_waves_over_a_bar",
        "SV_17_macdonald_rc_sub",
        "SV_18_macdonald_rv_sub",
        "SV_19_macdonald_rv_trans",
        "SV_20_dambreak_triangular_obstacle",
        "SV_21_dambreak_constriction",
        "SV_22_channel_constriction",
        "SV_23_cadam_lake_at_rest",
        "SV_24_cadam_steady_state",
        "SV_25_cadam_steady_state_friction",
        "SV_26_variable_friction",
        ]

    for case in caselist:
        rep = os.path.join(ARRAKIS_DIR, "examples", case)

        try:
            # Build input files if requested
            if ini:
                subprocess.run([sys.executable, "build_input_files.py"], cwd=rep, check=True)
            
            # Find and run test file
            test_files = glob.glob(os.path.join(rep, "*test.py"))
            if not test_files:
                print("Warning: No test file found in {}".format(case))
                continue
            
            test_file = os.path.basename(test_files[0])
            subprocess.run([sys.executable, test_file], cwd=rep, check=True)
            
        except subprocess.CalledProcessError:
            print("Check failed in study case: {}".format(case))
            return 1

    return 0            

def validate_sve(ini=False):
    """
    Run test located in examples
    """
    caselist = [\
        "SVE_00_bedload_dambreak",
        "SVE_01_bedload_bumpsub",
        "SVE_02_bedload_bumpcri",
        "SVE_03_bedload_bumptrans",
        "SVE_04_bedload_antidune",
        "SVE_05_bedload_lakeatrest",
        "SVE_06_bedload_weir",
        "SVE_07_bedload_stream",
        "SVE_08_bedload_constriction",
        "SVE_09_bedload_constriction_nacher",
        "SVE_10_bedload_cadam_lake_at_rest",
        "SVE_11_bedload_soni",
        "SVE_12_bedload_newton",
        "SVE_15_suspension_deposition",
        "SVE_16_suspension_vanrijn",
        "SVE_18_slope_effect",
        "SVE_22_bedload_sand_pile",
        ]

    for case in caselist:
        rep = os.path.join(ARRAKIS_DIR, "examples", case)

        try:
            # Build input files if requested
            if ini:
                subprocess.run([sys.executable, "build_input_files.py"], cwd=rep, check=True)
            
            # Find and run test file
            test_files = glob.glob(os.path.join(rep, "*test.py"))
            if not test_files:
                print("Warning: No test file found in {}".format(case))
                continue
            
            test_file = os.path.basename(test_files[0])
            subprocess.run([sys.executable, test_file], cwd=rep, check=True)
            
        except subprocess.CalledProcessError:
            print("Check failed in study case: {}".format(case))
            return 1

    return 0

def validate_st(ini=False):
    """
    Run test located in examples
    """
    caselist = [\
        "ST_00_bumpsub",
        "ST_01_flume_with_friction",
        "ST_02_macdonald_rv_sub",
        "ST_03_macdonald_rv_trans",
        "ST_04_tracer",
        "ST_05_tracer_multiple_rv",
        "ST_06_saint_matthieu",
        "STE_00_bedload_bumpsub",
        ]

    for case in caselist:
        rep = os.path.join(ARRAKIS_DIR, "examples", case)

        try:
            # Build input files if requested
            if ini:
                subprocess.run([sys.executable, "build_input_files.py"], cwd=rep, check=True)
            
            # Find and run test file
            test_files = glob.glob(os.path.join(rep, "*test.py"))
            if not test_files:
                print("Warning: No test file found in {}".format(case))
                continue
            
            test_file = os.path.basename(test_files[0])
            subprocess.run([sys.executable, test_file], cwd=rep, check=True)
            
        except subprocess.CalledProcessError:
            print("Check failed in study case: {}".format(case))
            return 1

    return 0

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-sv'  , default=False, action="store_true", help='run SV validation')
    parser.add_argument('-sve' , default=False, action="store_true", help='run SVE validation')
    parser.add_argument('-st'  , default=False, action="store_true", help='run ST validation')
    parser.add_argument('-ini' , default=False, action="store_true", help='build iniput files')
    ARGS = parser.parse_args()
    
    NO_ARGS = bool(ARGS.sv is False and ARGS.sve is False and ARGS.st is False)

    if NO_ARGS or ARGS.sv:
        print("~~> Validation of SV")
        output_sv = validate_sv(ini=ARGS.ini)
        if output_sv == 1:
            sys.exit("ARRAKIS SV VALIDATION FAILED")

    if NO_ARGS or ARGS.sve:
        print("~~> Validation of SVE")
        output_sve = validate_sve(ini=ARGS.ini)
        if output_sve == 1:
            sys.exit("ARRAKIS SVE VALIDATION FAILED")

    if NO_ARGS or ARGS.st:
        print("~~> Validation of ST")
        output_st = validate_st(ini=ARGS.ini)
        if output_st == 1:
            sys.exit("ARRAKIS ST VALIDATION FAILED")
