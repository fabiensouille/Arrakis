#!/usr/bin/env python
import numpy as np
from arrakis import *

ROOT_DIR = os.path.dirname(__file__)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nx', default=10000, type=int)
    args = parser.parse_args()

    # reference solution
    case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bumptrans.yml"))
    case_inputs['OUTPUT PARAMETERS']['LISTING FREQUENCY'] = 100000
    case_inputs['OUTPUT PARAMETERS']['OUTPUT FREQUENCY']  = 100000
    case_inputs['TIME PARAMETERS']['MAXIMAL NUMBER OF ITERATIONS'] = 10000000
    case_inputs['MESH PARAMETERS']['NUMBER OF NODES'] = args.nx
    case_inputs['NUMERICAL PARAMETERS']['CFL NUMBER'] = 0.9
    case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 4
    case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 1
    case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
    case_inputs['NUMERICAL PARAMETERS']['H FLUX LIMITOR'] = 1
    case_inputs['NUMERICAL PARAMETERS']['U FLUX LIMITOR'] = 1
    case_inputs['NUMERICAL PARAMETERS']['RECONSTRUCTION METHOD'] = 2
    case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
    run_arrakis(case_inputs)

    # save in .txt
    time, res, var = load_res(resname="RES", record=-1)
    
    # correction for exact flowrate
    for i in range(args.nx):
        q = 0.4
        res[i, 3] = q/res[i, 2]
    
    # save ref
    np.savetxt("REF/REFERENCE_SOLUTION_N{}.txt".format(args.nx),\
        np.c_[res[:,0], res[:,1], res[:,2], res[:,3]])
