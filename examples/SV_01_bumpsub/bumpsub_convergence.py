#!/usr/bin/env python
import os
import re
import numpy as np
from scipy import interpolate
import matplotlib.pylab as plt
from arrakis import *
from build_input_files import h as h_funct
from build_input_files import zb as zb_funct

ROOT_DIR = os.path.dirname(__file__)

def run_convergence(cas, lx, nx_list, schemes, scheme, torder, sorder, sol=None, ref=None):
    """
    Run convergence test
    """
    if sol is None and ref is None:
        raise ValueError("Provide either ref of sol")

    n_sc = len(schemes)
    n_mesh = len(nx_list)
    dx_list = [lx/i for i in nx_list]

    # error lists:
    error_L1_h = np.empty((n_sc, n_mesh))
    error_L1_u = np.empty((n_sc, n_mesh))
    error_L1_q = np.empty((n_sc, n_mesh))
    error_L2_h = np.empty((n_sc, n_mesh))
    error_L2_u = np.empty((n_sc, n_mesh))
    error_L2_q = np.empty((n_sc, n_mesh))

    # loop on schemes:
    for j, sc in enumerate(schemes):

        space_scheme = scheme[j]
        time_scheme = torder[j]
        order = sorder[j]

        # loop on meshes:
        for m, nx in enumerate(nx_list):

            # read case.yml
            case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, cas))

            # build input files
            build_geo_file(nx=nx, xa=0., xb=10., 
                           geo_file="MESH_TMP.geo", zb_funct=zb_funct)
            build_ini_file(geo_file="MESH_TMP.geo", 
                           ini_file="INIC_TMP.ini", h_funct=h_funct, q=0.2)

            # modify options
            case_inputs['TIME PARAMETERS']['MAXIMAL NUMBER OF ITERATIONS'] = 10000000
            case_inputs['TIME PARAMETERS']['DURATION'] = 1000.
            case_inputs['MESH PARAMETERS']['MESH FILE NAME'] = "MESH_TMP"
            case_inputs['INITIAL CONDITION']['INITIAL FILE NAME'] = "INIC_TMP"
            case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = space_scheme
            case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = time_scheme
            case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = order
            case_inputs['NUMERICAL PARAMETERS']['H FLUX LIMITOR'] = 1
            case_inputs['NUMERICAL PARAMETERS']['U FLUX LIMITOR'] = 1
            case_inputs['NUMERICAL PARAMETERS']['RECONSTRUCTION METHOD'] = 2
            case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1

            # run arrakis
            run_arrakis(case_inputs)

            # get res:
            time, res, var = load_res(resname="RES", record=-1)
            x = res[:,0]
            zb= res[:,1]
            h = res[:,2]
            u = res[:,3]

            # ref:
            if sol is not None:
                # compute ref:
                sol()
                # interpolate analyticsol on x
                hinterp = interpolate.interp1d(sol.x, sol.H)
                uinterp = interpolate.interp1d(sol.x, sol.U)
                
            elif ref is not None:
                # interpolate analyticsol on x
                hinterp = interpolate.interp1d(ref[:, 0], ref[:, 2]) #H
                uinterp = interpolate.interp1d(ref[:, 0], ref[:, 3]) #U
            
            href = hinterp(x)
            uref = uinterp(x)
            qref = href*uref

            # compute L1,L2 errors:
            error_L1_h[j, m] = error_L1(h, href)
            error_L1_u[j, m] = error_L1(u, uref)
            error_L1_q[j, m] = error_L1(h*u, href*uref)
            error_L2_h[j, m] = error_L2(h, href)
            error_L2_u[j, m] = error_L2(u, uref)
            error_L2_q[j, m] = error_L2(h*u, href*uref)

            # plot res vs ref
            #plot_res(time, res, var, "FIG/plot.png", ref=ref, show=True)

    # cleaning:
    os.system("rm RESU/*.dat")
    os.system("rm MESH_TMP.geo")
    os.system("rm INIC_TMP.ini")

    # save error files
    np.savetxt("RESU/error_L1_h.dat", error_L1_h)
    np.savetxt("RESU/error_L1_u.dat", error_L1_u)
    np.savetxt("RESU/error_L1_q.dat", error_L1_q)
    np.savetxt("RESU/error_L2_h.dat", error_L2_h)
    np.savetxt("RESU/error_L2_u.dat", error_L2_u)
    np.savetxt("RESU/error_L2_q.dat", error_L2_q)


def plot(lx, nx_list, schemes, scheme, torder, sorder):
    """
    Plot convergence results
    """
    n_sc = len(schemes)
    n_mesh = len(nx_list)
    dx_list = [lx/i for i in nx_list]
    xlim=[0.01, 1]

    # read error files
    error_L2_h = np.loadtxt("RESU/error_L2_h.dat")
    error_L2_u = np.loadtxt("RESU/error_L2_u.dat")
    error_L2_q = np.loadtxt("RESU/error_L2_q.dat")
    n_sc = np.shape(error_L2_h)[0]

    # plot errors:
    plot_convergence(
        [dx_list for i in range(n_sc)],
        [error_L2_h[i, :] for i in range(n_sc)],
        figsize=(5, 4),
        label=schemes,
        xlim=xlim,
        xlabel="$\Delta x$",
        ylabel="$E_h$ in L2 norm",
        figname="FIG/L2_error_h.png",
        dpi=300)

    plot_convergence(
        [dx_list for i in range(n_sc)],
        [error_L2_q[i, :] for i in range(n_sc)],
        figsize=(5, 4),
        label=schemes,
        xlim=xlim,
        xlabel="$\Delta x$",
        ylabel="$E_q$ in L2 norm",
        figname="FIG/L2_error_q.png",
        dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', default=False, action="store_true", help='run arrakis')
    parser.add_argument('-p', '--post', default=False, action="store_true", help='post processing')
    args = parser.parse_args()

    if args.run==False and args.post==False:
        no_args = True
    else:
        no_args = False
    
    # convergence options
    nx = [25, 50, 100, 200, 400, 800]
    names = ['Roe-Euler-1', 'Roe-Heun-2', 'HLL-Euler-1', 'HLL-Heun-2', 'Kin-Euler-1', 'Kin-Heun-2']
    scheme = [1, 1, 2, 2, 4, 4]
    torder = [1, 2, 1, 2, 1, 2]
    sorder = [1, 2, 1, 2, 1, 2]
   
    # reference solution
    sol = BumpAnalyticSol(flow='sub', Q=0.2, hl=0.4, bottom_function='exponential')
    #ref = np.loadtxt("REF/REFERENCE_SOLUTION_N10000.txt")
    
    # run 
    if no_args or args.run:
        run_convergence("bumpsub.yml", 10., nx, names, scheme, torder, sorder, sol=sol)
        #run_convergence("bumpsub.yml", 10., nx, names, scheme, torder, sorder, ref=ref)

    # post
    if args.post:
        plot(10., nx, names, scheme, torder, sorder)
