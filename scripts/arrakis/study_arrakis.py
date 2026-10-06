#!/usr/bin/env python3
"""
Study class
"""
import os
import sys
import glob
import shutil
import argparse
import numpy as np
from arrakis import read_yml_input_file, run_arrakis

class Study():
    """
    Study class

    Principle: can be used to run multiple cases for vnv or studies
    """
    def __init__(self):
        self.name = "study"
        self.cases = []
        self.global_dir = ""
        self.local_dir = ""

    def set_dir(self, file=__file__):
        """ set directories path """
        self.global_dir = os.path.dirname(os.path.dirname(os.path.abspath(file)))
        self.local_dir = os.path.dirname(file)

    def pre(self):
        """ Pre """
        self.set_dir()
        case_inputs = read_yml_input_file(os.path.join(self.global_dir, "scripts", "arrakis.yml"))
        self.cases.append(case_inputs)

    def run(self, omp_num_threads=4):
        """ Run """
        rep = os.path.join(self.local_dir, "RESU")
        # Remove all .dat files in RESU directory
        for file in glob.glob(os.path.join(rep, "*.dat")):
            if os.path.isfile(file):
                os.remove(file)
        for case in self.cases:
            run_arrakis(case, omp_num_threads=omp_num_threads)

    def check(self):
        """ Check """
        for case in self.cases:
            resf = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_fin = os.path.join(self.local_dir, 'RESU', resf + 'fin.dat')
            file_ref = os.path.join(self.local_dir, 'REF', resf + 'ref_{}.txt'.format(self.name))
            res = np.loadtxt(file_fin, delimiter=',', skiprows=2)
            ref = np.loadtxt(file_ref, delimiter=',', skiprows=2)
            try:
                np.testing.assert_allclose(res, ref, rtol=1e-7)
            except AssertionError:
                sys.exit("Check failed in study: {}, case: {}".format(self.name, resf))

    def reset_ref(self):
        """ Set reference file """
        os.makedirs("REF", exist_ok=True)

        for case in self.cases:
            res = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_fin = os.path.join(self.local_dir, 'RESU', res + 'fin.dat')
            file_ref = os.path.join(self.local_dir, 'REF', res + 'ref_{}.txt'.format(self.name))
            shutil.copy2(file_fin, file_ref)

    def post(self):
        """ Post """

    def execute(self, args):
        """
        Execute study based on input args
        """
        self.pre()

        no_args = bool(args.run is False and \
                       args.post is False and \
                       args.check is False and \
                       args.reset_ref is False)

        if no_args:
            self.run(omp_num_threads=args.ncsize)
            self.check()
            self.post()
        else:
            if args.run:
                self.run(omp_num_threads=args.ncsize)
            if args.check:
                self.check()
            if args.post:
                self.post()
            if args.reset_ref:
                self.reset_ref()

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', default=False, action="store_true", help='run arrakis')
    parser.add_argument('-p', '--post', default=False, action="store_true", help='post processing')
    parser.add_argument('-c', '--check', default=False, action="store_true", help='check results')
    parser.add_argument('--reset-ref', default=False, action="store_true", help='WARNING: only use in extreme necessity')
    parser.add_argument('-nc', '--ncsize', type=int, default=4, help='number of OpenMP threads (default: 4)')
    ARGS = parser.parse_args()

    # test
    study = Study()
    study.execute(ARGS)
