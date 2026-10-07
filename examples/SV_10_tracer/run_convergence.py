import os
import numpy as np
from scipy import interpolate
import matplotlib.pylab as plt
from arrakis import *
from build_input_files import tracers

ROOT_DIR = os.path.dirname(__file__)

# convergence options: 'SPACE' or 'TIME'
SCHEMES = 'SPACE'
CONVERGENCE = 'SPACE'

# parameters:
PLOT  = False
CLEANING = True

def main():
    """
    Main function
    """
    # schemes:
    if SCHEMES=='SPACE':
        schemes = ['KI1', 'KI2']
        schemes_id = [5, 5]
        sorder = [1, 2]
        colors = ["b", "g"]
        n_sc = len(schemes)

    elif SCHEMES=='TIME':
        schemes = ['EUL', 'HEU']
        schemes_id = [1, 2]
        colors = ["b", "g"]
        n_sc = len(schemes)

    else:
        raise ValueError("Unknown convergence type")

    # discretisation:
    if CONVERGENCE=='SPACE':
        nx_list = [51, 201, 401, 801]
        dx_list = [10./i for i in nx_list]
        n_mesh = len(nx_list)
        dt_list = [0.9 for i in range(n_mesh)]
        xlim=[0.01, 1]
        step_list = dx_list

    elif CONVERGENCE=='TIME':
        dt_list = [-1.e-2, -0.5e-2, -1.e-3, -0.5e-3]
        n_mesh = len(dt_list)
        nx_list = [51 for i in range(n_mesh)]
        xlim=None
        step_list = [abs(dt_list[i]) for i in range(n_mesh)]

    else:
        raise ValueError("Unknown convergence type")

    # error lists:
    error_L1_T = np.empty((n_sc, n_mesh))
    error_L2_T = np.empty((n_sc, n_mesh))

    # loop on schemes:
    for j, sc in enumerate(schemes):

        if SCHEMES=='SPACE':
            time_scheme = 2
            space_scheme = schemes_id[j]
            order = sorder[j]
        elif SCHEMES=='TIME':
            time_scheme = schemes_id[j]
            space_scheme = 5
            order = 2

        # loop on meshes:
        for m, nx in enumerate(nx_list):

            # read case.yml
            case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "tracer_canal.yml"))

            # build input files:
            build_geo_file(nx=nx, xa=0., xb=10., geo_file="MESH_TMP.geo")
            build_ini_file(geo_file="MESH_TMP.geo",
                           ini_file="INIC_TMP.ini", h=1., u=0.1,
                           ntrac=1, t_funct=tracers)

            # modify options
            case_inputs['MESH PARAMETERS']['MESH FILE NAME'] = "MESH_TMP"
            case_inputs['INITIAL CONDITION']['INITIAL FILE NAME'] = "INIC_TMP"
            case_inputs['TIME PARAMETERS']['CFL NUMBER'] = dt_list[m]
            case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = space_scheme
            case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = time_scheme
            case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = order

            # run arrakis
            run_arrakis(case_inputs)
       
            # get res:
            time, res, _ = load_res(resname="RES", record=-1)
            rest = np.loadtxt("RESU/RES_TRAC0001_fin.dat", delimiter=',', skiprows=1)

            x = rest[:,0]
            T = rest[:,1]

            # compute ref:
            sol = TraceurAdvection(u0=0.1, length=10.)
            sol(time)
            ref = [sol.x, sol.T]

            # interpolate analyticsol on x
            Tinterp = interpolate.interp1d(sol.x, sol.T)
            Tref = Tinterp(x)

            # compute L1 errors:
            error_L1_T[j, m] = error_L1(T, Tref)

            # compute L2 errors:
            error_L2_T[j, m] = error_L2(T, Tref)

            # plot res vs ref
            if PLOT:
                plot_res(res, "FIG/plot.png", ref=ref, show=True)

            # cleaning:
            if CLEANING:
                os.system("rm RESU/*")
                os.system("rm *_TMP*")

    # plot errors:
    figname = "FIG/L2_error_T_{}schemes_vs_{}step.png".format(
        SCHEMES.lower(), CONVERGENCE.lower())

    plot_convergence(
        [step_list for i in range(n_sc)],
        [error_L2_T[i, :] for i in range(n_sc)],
        label=schemes,
        xlim=xlim,
        figname=figname,
        dpi=300)

if __name__ == '__main__':
    main()
