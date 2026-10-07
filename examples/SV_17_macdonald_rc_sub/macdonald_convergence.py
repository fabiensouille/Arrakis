#!/usr/bin/env python
import os
import re
import numpy as np
from scipy import interpolate
import matplotlib.pylab as plt
from arrakis import *
from build_input_files import bottom, depth, velocity

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

    # loop on meshes:
    for m, nx in enumerate(nx_list):

        # create mesh and ini files
        build_geo_file(
            nx=nx, xa=0., xb=1000.,
            geo_file="MESH_TMP.geo", 
            zb_funct=bottom,
            variable_width=True, 
            Bx=10.)
        build_ini_file(
            geo_file="MESH_TMP.geo",
            ini_file="INIC_TMP.ini",
            h_funct=depth,
            u_funct=velocity)

        # loop on schemes:
        for j, sc in enumerate(schemes):

            space_scheme = scheme[j]
            time_scheme = torder[j]
            order = sorder[j]

            # read case.yml
            case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, cas))

            # modify options
            case_inputs['MESH PARAMETERS']['MESH FILE NAME'] = 'MESH_TMP'
            case_inputs['INITIAL CONDITION']['INITIAL FILE NAME'] = 'INIC_TMP'
            case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = space_scheme
            case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = time_scheme
            case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = order

            # run arrakis
            run_arrakis(case_inputs)

            # get res:
            time, res, var = load_res(resname="RES", record=-1)
            x = res[:,0]
            zb= res[:,1]
            h = res[:,2]
            u = res[:,3]
            L = res[:,4]
            q = res[:,5]

            # ref:
            if sol is not None:
                # interpolate analyticsol on x
                hinterp = interpolate.interp1d(sol.x, sol.H)
                uinterp = interpolate.interp1d(sol.x, sol.U)
                qinterp = interpolate.interp1d(sol.x, sol.Q)
            elif ref is not None:
                # interpolate analyticsol on x
                hinterp = interpolate.interp1d(ref[:, 0], ref[:, 2]) #H
                uinterp = interpolate.interp1d(ref[:, 0], ref[:, 3]) #U
                qinterp = interpolate.interp1d(ref[:, 0], ref[:, 5]) #Q

            href = hinterp(x)
            uref = uinterp(x)
            qref = qinterp(x)

            # compute L1,L2 errors:
            error_L1_h[j, m] = error_L1(h, href)
            error_L1_u[j, m] = error_L1(u, uref)
            error_L1_q[j, m] = error_L1(q, qref)
            error_L2_h[j, m] = error_L2(h, href)
            error_L2_u[j, m] = error_L2(u, uref)
            error_L2_q[j, m] = error_L2(q, qref)

            # plot res vs ref
            #plot_res(time, res, var, "FIG/plot.png", ref=ref, show=True)

    # cleaning:
    os.system("rm RESU/*.dat")
    os.system("rm *_TMP*")

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
    #sol = McDonnaldRCSubAnalyticSol(xa=0., xb=1000.)
    #sol()
    ref = np.loadtxt("REF/ANALYTICAL_SOLUTION_N20000.txt")

    # run 
    if no_args or args.run:
        #run_convergence("macdonald_rc_sub.yml", 10., nx, names, scheme, torder, sorder, sol=sol)
        run_convergence("macdonald_rc_sub.yml", 10., nx, names, scheme, torder, sorder, ref=ref)

    # post
    if args.post:
        plot(10., nx, names, scheme, torder, sorder)
