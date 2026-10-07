#!/usr/bin/env python3
"""
Basic plots of results

How to execute: plot.py <my_case>.yml
Note : you should have run your case first with: arrakis.py <my_case>.yml
"""
import os
import argparse
from post import load_res, plot_res
from arrakis import read_yml_input_file

if __name__ == '__main__':

    # parse options
    parser = argparse.ArgumentParser(description='Plot arrakis')
    parser.add_argument('case', type=str, nargs='+', help='input arrakis data file (.yml)')
    parser.add_argument('-b', '--postb', default=False,
                        action="store_true", help='post processing')
    parser.add_argument('-i', '--irecord', default=-1, type=int, help='record to plot')
    args = parser.parse_args()

    # read yaml input file
    case_dict = read_yml_input_file(args.case[0])
    resname = case_dict['OUTPUT PARAMETERS']['RESULT FILE']

    # post processing for SVE
    if args.postb:
        time, res, var = load_res(resname=resname, record=args.irecord)
        plot_res(time, res, var,\
                 os.path.join("FIG", "plotb.png"), show=True, plot_sediments=True)

    # post processing for SV
    else:
        time, res, var = load_res(resname=resname, record=args.irecord)
        plot_res(time, res, var,\
                 os.path.join("FIG", "plot.png"), show=True)
