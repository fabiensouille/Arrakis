#!/usr/bin/env python
from arrakis import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nx', default=10000, type=int)
    args = parser.parse_args()

    # analytical solution
    sol = McDonnaldRVTransAnalyticSol(xa=0., xb=1000., N=args.nx)
    sol()

    # save in .txt
    sol.savetxt("REF/ANALYTICAL_SOLUTION_N{}.txt".format(args.nx))
