import numpy as np
from .mesh import Mesh1D, refine_mesh
from .utils import smooth_data
from scipy.interpolate import interp1d

def build_geo_file(geo_file="MESH.geo",
                   mesh=None, nx=100, xa=0., xb=10., density=None,
                   erodible_bed=False, variable_width=False,
                   zb= 0., zb_funct=None, zb_vect=None, x_vect=None,
                   Bx= 1., Bx_funct=None, Bx_vect=None,
                   zf=-1., zf_funct=None, zf_vect=None,
                   verbose=False):
    """ 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Build geometry file : contains only internal nodes x[1,nx]

        xA    x1               xi              xn   xB   
    o----|----o--  ...  --|----o----|--  ...  --o----|----o
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param geo_file (str): geometry file (default: "mesh.geo")
    @param nx (int): number of cells (default: 100)
    @param xa (int): abscissa of left boundary interface (default: 0.)
    @param xb (int): abscissa of right boundary interface (default: 10.)
    @param density[0:nx+1] (float): array containing density (default: None)
                                    used to define cell size
    @param erodible_bed (bool): write zf in mesh file (default: False)
    @param variable_width (bool): write width Bx in mesh file (default: False)
    @param zb (float): constant bottom (default: 0.)
    @param zb_funct (function): bottom function (default: None)
    @param zb_vect (array): bottom vector (default: None)
    @param zf (float): constant hard bottom (default: -1.)
    @param zf_funct (function): hard bottom function (default: None)
    @param zf_vect (array): hard bottom vector (default: None)
    @param Bx (float): constant width (default: 1.)
    @param Bx_funct (function): width function (default: None)
    @param Bx_vect (array): width vector (default: None)
    @param x_vect (array): abscissa vector for zb, zf, Bx (default: None)
    @param verbose (bool): verbose mode (default: False)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # build 1d mesh
    if mesh is None:
        mesh1d = Mesh1D(nx=nx, xa=xa, xb=xb, density=density, verbose=verbose)
    else:
        mesh1d = mesh
    x = mesh1d.x[1:-1]

    # check vector intputs and create interpolation functions 
    if zb_vect is not None and len(zb_vect)!=nx:  
        if x_vect is None:
            raise ValueError("x_vect required if zb_vect given and len(zb_vect)!=nx")
        fzb = interp1d(x_vect, zb_vect, kind='linear', fill_value="extrapolate")
    else:
        fzb = None
    if zf_vect is not None and len(zf_vect)!=nx:
        if x_vect is None:
            raise ValueError("x_vect required if zf_vect given and len(zf_vect)!=nx")
        fzf = interp1d(x_vect, zf_vect, kind='linear', fill_value="extrapolate")
    else:
        fzf = None
    if Bx_vect is not None and len(Bx_vect)!=nx:
        if x_vect is None:
            raise ValueError("x_vect required if Bx_vect given and len(Bx_vect)!=nx")
        fBx = interp1d(x_vect, Bx_vect, kind='linear', fill_value="extrapolate")
    else:
        fBx = None

    # initialize arrays
    zb_arr = np.empty(nx, dtype='d')
    zf_arr = np.empty(nx, dtype='d')
    Bx_arr = np.empty(nx, dtype='d')

    # init bottom and width - loop on mesh nodes
    for i in range(nx):

        # create sediment bed
        if zb_funct is not None:
            zb_arr[i] = zb_funct(x[i])
        elif zb_vect is not None:
            if fzb is not None:
                zb_arr[i] = fzb(x[i])
            else:
                zb_arr[i] = zb_vect[i]
        else:
            zb_arr[i] = zb

        # create non erodible bed
        if erodible_bed:
            if zf_funct is not None:
                zf_arr[i] = zf_funct(x[i])
            elif zf_vect is not None:
                if fzf is not None:
                    zf_arr[i] = fzf(x[i])
                else:
                    zf_arr[i] = zf_vect[i]
            else:
                zf_arr[i] = zf

        # create width
        if variable_width:
            if Bx_funct is not None:
                Bx_arr[i] = Bx_funct(x[i])
            elif Bx_vect is not None:
                if fBx is not None:
                    Bx_arr[i] = fBx(x[i])
                else:
                    Bx_arr[i] = Bx_vect[i]
            else:
                Bx_arr[i] = Bx

    # write .geo file
    if variable_width:
        write_geo_file(geo_file, x, zb_arr, ZF=None, L=Bx_arr)
        if erodible_bed:
            write_geo_file(geo_file, x, zb_arr, ZF=zf_arr, L=Bx_arr)
    else:
        write_geo_file(geo_file, x, zb_arr, ZF=None, L=None)
        if erodible_bed:
            write_geo_file(geo_file, x, zb_arr, ZF=zf_arr, L=None)

    # verbose            
    if verbose:
        print("zb = :", zb_arr)
        if erodible_bed:
            print("zf = :", zf_arr)
        if variable_width:
            print("Bx = :", Bx_arr)

    return 0

def write_geo_file(geo_file, X, ZB, ZF=None, L=None):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Write Arrakis geometry file
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param geo_file (str): geometry file (default: "mesh.geo")
    @param X (np.array[nx], float): x abscissa of profiles
    @param ZB (np.array[nx], float): bed elevation
    @param ZF (np.array[nx], float): non erodible bed elevation
    @param L (np.array[nx], float): width
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # number of points
    nx = len(X)
    assert len(ZB)==nx

    # write .geo file
    f = open(geo_file, "w")
    f.write("# Mesh file, Nx={} \n".format(nx))

    if ZF is None and L is None:
        f.write("X,ZB \n")
        for i in range(nx):
            f.write("{:.8f} {:}\n".format(X[i], ZB[i]))

    elif ZF is None and L is not None:
        assert len(L)==nx
        f.write("X,ZB,L \n")
        for i in range(nx):
            f.write("{:.8f} {:} {:}\n".format(X[i], ZB[i], L[i]))

    elif ZF is not None and L is None:
        assert len(ZF)==nx
        f.write("X,ZB,ZF \n")
        for i in range(nx):
            f.write("{:.8f} {:} {:}\n".format(X[i], ZB[i], ZF[i]))

    elif ZF is not None and L is not None:
        f.write("X,ZB,ZF,L \n")
        assert len(ZF)==nx
        assert len(L)==nx
        for i in range(nx):
            f.write("{:.8f} {:} {:} {:}\n".format(X[i], ZB[i], ZF[i], L[i]))
    f.close()

    return 0

def read_geo_file(geo_file):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read Arrakis geometry file
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param geo_file (str): geometry file (default: "mesh.geo")
    """
    # open file
    f = open(geo_file, "r")

    # read number of nodes
    line = f.readline().split(",")
    nx = int(line[1].split("=")[1])

    # read header and check content
    line = f.readline().split(",")
    line = [item.strip() for item in line]
    if "X" not in line or "ZB" not in line:
        raise ValueError("geo file must contain X and ZB")
    f.close()

    # read data
    data = np.loadtxt(geo_file, delimiter=" ", skiprows=2)
    X = data[:, line.index("X")]
    ZB = data[:, line.index("ZB")]
    if "ZF" in line:
        ZF = data[:, line.index("ZF")]
    else:
        ZF = None
    if "L" in line:
        L = data[:, line.index("L")]
    else:
        L = None

    return X, ZB, ZF, L

def read_mascaret_geo(file, verbose=True, no_data_value=-9999.):
    """ 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read Mascaret .geo file
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param file (str): Mascaret geo file

    @return x (np.array[nprof], float): x abscissa of profiles
    @return nprof (int): number of profiles
    @return ny (np.array[nprof], int): number of points in each profiles
    @return profiles (np.array[nprof, np.max(ny), 2], float): 
        Y, Zb data for each profile
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # First read : count number and size of profiles
    # ~~~~~~~~~~
    nprof = 0 # Number of profiles in geoC
    ny = []   # Number of points along y axis in each profiles
    cny = 0   # Counter for ny
    f = open(file, "r")
    for line in f:
        line_content = line.split()
        if (line_content[0]=='Profil' or line_content[0]=='PROFIL'):
            if nprof!=0:
               ny.append(cny)
            nprof += 1
            cny = 0
        else:
            cny += 1
    ny.append(cny)
    ny = np.asarray(ny)
    assert len(ny)==nprof
    f.close()

    # Second read : read file content and append profiles
    # ~~~~~~~~~~~
    x = np.zeros((nprof), dtype='d')
    profiles = no_data_value*np.ones((nprof, np.max(ny), 2), dtype='d')
    iprof = -1 # Index of profile
    cny = 0    # Counter for ny
    f = open(file, "r")
    for line in f:
        line_content = line.split()
        if (line_content[0]=='Profil' or line_content[0]=='PROFIL'):
            iprof += 1
            x[iprof] = line_content[3]
            cny = 0
        else:
            profiles[iprof, cny, 0] = float(line_content[0])
            profiles[iprof, cny, 1] = float(line_content[1])
            cny += 1
    f.close()

    if verbose:
        print("nprof =", nprof)
        print("ny =", ny)
        print("profiles =", profiles)
        print("x =", x)

    return x, nprof, ny, profiles

def read_courlis_geoC(file, verbose=True, no_data_value=-9999.):
    """ 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read Courlis .geoC file 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param file (str): courlis geoC file

    @return x (np.array[nprof], float): x abscissa of profiles
    @return nprof (int): number of profiles
    @return ny (np.array[nprof], int): number of points in each profiles
    @return profiles (np.array[nprof, np.max(ny), 3], float): 
        Y, Zb, Zf data for each profile
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # First read : count number and size of profiles
    # ~~~~~~~~~~
    nprof = 0 # Number of profiles in geoC
    ny = []   # Number of points along y axis in each profiles
    cny = 0   # Counter for ny
    f = open(file, "r")
    for line in f:
        line_content = line.split()
        if (line_content[0]=='Profil' or line_content[0]=='PROFIL'):
            if nprof!=0:
               ny.append(cny)
            nprof += 1
            cny = 0
        else:
            cny += 1
    ny.append(cny)
    ny = np.asarray(ny)
    assert len(ny)==nprof
    f.close()

    # Second read : read file content and append profiles
    # ~~~~~~~~~~~
    x = np.zeros((nprof), dtype='d')
    profiles = no_data_value*np.ones((nprof, np.max(ny), 3), dtype='d')
    iprof = -1 # Index of profile
    cny = 0    # Counter for ny
    f = open(file, "r")
    for line in f:
        line_content = line.split()
        if (line_content[0]=='Profil' or line_content[0]=='PROFIL'):
            iprof += 1
            x[iprof] = line_content[3]
            cny = 0
        else:
            profiles[iprof, cny, 0] = float(line_content[0])
            profiles[iprof, cny, 1] = float(line_content[1])
            profiles[iprof, cny, 2] = float(line_content[2])
            cny += 1
    f.close()

    if verbose:
        print("nprof =", nprof)
        print("ny =", ny)
        print("profiles =", profiles)
        print("x =", x)

    return x, nprof, ny, profiles

def convert_profile(y, zb_in, zf_in=None, zs_in=None, h_in=None, 
                    option=3, option_z=2, eps=1.e-6, 
                    verbose=False, plot_profile=False):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Convert arbitrary profile into rectangular profile
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param y (np.array[ny], float): profile abscissa 
    @param zb_in (np.array[ny], float): bed elevation
    @param zf_in (np.array[ny], float): non erodible bed elevation
    @param zs_in (float): free surface elevation
    @param option (int): option for convertion into rectangular cross-section
        option=1 : zb=Zb, zs=Zs and B=A/(Zs-Zb) 
        option=2 : zb=Zb, B=B(Mascaret) computed where Z_s-Z_b>0
        option=3 : zb=Zb, zs=Zb+A/B, B such that Rh=Rh(Mascaret) (default)
    @param option_z (int): option for vertical position of rectangular profile
        option_z=1 : zb=min(Zb(y))
        option_z=2 : zb=mean(Zb(y)) (default)
        option_z=3 : zs=Zs
    @param eps (float): epsilon to define wet area (default: 1.e-8)
    @param verbose (bool): verbose mode (default: False)
    @param plot_profile (bool): plot profile (default: False)

    @return B (float): profile width
    @return zb_out (float): profile bed elevation
    @return zf_out (float): profile non erodible bed elevation
    @return zs_out (float): free surface elevation
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # Initialize parameters
    # ~~~~~~~~~~~~~~~~~~~~~
    if zs_in is not None and h_in is not None:
        raise ValueError("set zs_in OR h_in")

    # if no surface data: get max of bottom elevation
    if zs_in is None and h_in is None:
        zs_out = np.max(zb_in)
    else:
        if h_in is not None:
            zs_out = np.min(zb_in) + h_in
        else:
            zs_out = zs_in

    zs0 = zs_out
    zb_out = np.min(zb_in)
    mean_width = np.mean(y)

    A = 0. # wet area
    P = 0. # wetted perimeter
    As = 0. # sediment area
    mean_zby = 0. # mean bed elevation over wet area
    leng_zby = 0

    y0_surf = -999. # first wet point of a pond
    y1_surf = -999. # last wet point of a pond
    is_wet = False
    wet_width = 0.

    if zf_in is not None:
        zf_out = np.min(zf_in)
        sediments = True
        y0_bed = -999.
        y1_bed = -999.
        is_mobile = False
        erodible_width =0.
    else:
        zf_out = zb_out
        sediments = False

    ny = len(y) # number of points in profile

    # loop on profile point along y
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # -> compute wet width and area A
    # -> compute erodible width and area As
    # -> compute mean zb over wet area
    for j in range(ny-1):
        zb = zb_in[j]
        condition = (zs_out-zb_in[j+1]>eps) or (zs_out-zb_in[j]>eps)

        # increment mean zb
        if condition:
            mean_zby += zb*(y[j+1]-y[j])
            leng_zby += y[j+1]-y[j]

        # Compute A_in
        if condition:
            A += (y[j+1]-y[j])*max(0., (0.5*(max(zs_out, zb_in[j+1]) \
                                            +max(zs_out, zb_in[j])) \
                                       -0.5*(zb_in[j+1]+zb_in[j])))
        # Compute P_in
        if condition:
            P += np.sqrt((y[j+1]-y[j])**2 + (zb_in[j+1]-zb_in[j])**2)

        # Compute As_in
        if sediments==True:
            As += (y[j+1]-y[j])*max(0.,(0.5*(zb_in[j+1]+zb_in[j]) \
                                       -0.5*(zf_in[j+1]+zf_in[j])))

    # finalize mean zb
    if leng_zby>0:
        mean_zby /= leng_zby

    for j in range(ny):
        zb = zb_in[j]
        zs = max(zs_out, zb)

        # Compute B_in
        if zs-zb > eps:
            if is_wet == False:
                is_wet = True
                y1_surf = -999.
                if y0_surf==-999.:
                    y0_surf = y[j]

        if zs-zb < eps or j==ny-1:
            if is_wet == True:
                is_wet = False
                if y1_surf==-999.:
                    y1_surf = y[j]
                wet_width += max(0., y1_surf-y0_surf)
                y0_surf = -999.

        if sediments==True:
            zf = zf_in[j]

            # Compute Bs_in
            if zb-zf > eps:
                if is_mobile == False:
                    is_mobile = True
                    y1_bed = -999.
                    if y0_bed==-999.:
                        y0_bed = y[j]

            if zb-zf < eps:
                if is_mobile == True:
                    is_mobile = False
                    if y1_bed==-999.:
                        y1_bed = y[j]
                    erodible_width += max(0., y1_bed - y0_bed)
                    y0_bed = -999.
    
    if verbose:
        print("Input cross section properties: ")
        print("A  =", A)
        print("P  =", P)
        print("B  =", wet_width)
        print("Zb =", min(zb_in))
        if sediments==True:
            print("As =", As)
            print("Bs =", erodible_width)

    # Compute rectangular profile
    # ***************************
    # B = A/(zs-zb)
    if option==1:
        if wet_width > eps:
            B = A/max(eps, zs_out-zb_out)
        else:
            B = mean_width

    # B = B
    elif option==2:
        if wet_width > eps:
            B = wet_width
        else:
            B = mean_width

    # B such that Rh=Rh
    elif option==3:
        Delta = P**2 - 8.*A 
        if wet_width > eps:
            if Delta>=0.:
                B = 0.5*(P + np.sqrt(Delta))
            else:
                B = A/max(eps, zs_out-zb_out)
        else:
            B = mean_width

    # Compute zs and zf to conserve A and As
    # **************************************
    zs_out = zb_out + A/max(eps, B)

    if sediments==True:
        zf_out = zb_out - As/max(eps, B)

    # Vertical translation of profile
    # *******************************
    if option_z==1:
        dz = 0.
    elif option_z==2:
        dz = mean_zby - zb_out
    elif option_z==3:
        dz = zs0 - zs_out

    zf_out += dz
    zb_out += dz
    zs_out += dz

    # debug print
    if verbose:
        print("Rectangular cross section properties: ")
        print("A  =", B*(zs_out-zb_out))
        print("P  =", B + 2.*(zs_out-zb_out))
        print("B  =", B)
        print("zb =", zb_out)
        print("zs =", zs_out)
        print("z shift of zb =", zb_out - min(zb_in))
        print("z shift of zs =", zs_out - zs_in)

    # plot 
    if plot_profile:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(9., 4.5))
        y_shift = 0.5*(max(y)-min(y)) - 0.5*B
        ax.plot(y, [max(zs_in, zb_in[j]) for j in range(len(y))], 
                marker='s', label='$Z_s$ (In)', c='b', markersize=3)
        ax.plot(y, zb_in, marker='s', label='$Z_b$ (In)', c='r', markersize=3)
        ax.fill_between(y, zb_in, zs_in, where=zs_in-zb_in>0, color='steelblue', alpha=0.25)
        ax.plot(y_shift+[0., B], [zs_out, zs_out], marker='o', label='$z_s$ (out)', 
                c='b', ls='--', markersize=3)
        ax.plot(y_shift+[0., B], [zb_out, zb_out], marker='o', label='$z_b$ (out)', 
                c='r', ls='--', markersize=3)
        ax.fill_between(y_shift+[0., B], zb_out, zs_out, color='grey', alpha=0.25)
        plt.legend()
        ax.grid()
        plt.show()

    return B, zb_out, zf_out, zs_out

def convert_geo_mascaret_to_arrakis(
        mascaret_file,
        arrakis_file="MESH.geo",
        xs_in=None,
        zs_in=None,
        h_in=None,
        option=3,
        option_z=2,
        verbose=False,
        plot_profile=False):
    """ 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Convert Mascaret .geo file into Arrakis .geo file
    
    Remark: convertion to rectangular profile requires water elevation data:
      - zs_in : array containing water surface elevation at xs_in 
      or 
      - h_in : same constant water depth for all profiles
      - if None of them is given, zs is set to max(zb)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param mascaret_file (str): Mascaret geo file
    @param arrakis_file (str): Arrakis geo file (default: "MESH.geo")
    @param xs_in (np.array[nx], float): x abscissa of zs_in
    @param zs_in (np.array[nx], float): free surface elevation
    @param h_in (float): water depth
    @param option (int): option for convertion into rectangular cross-section
        option=1 : zb=Zb, zs=Zs and B=A/(Zs-Zb)
        option=2 : zb=Zb, B=B(Mascaret) computed where Z_s-Z_b>0
        option=3 : zb=Zb, zs=Zb+A/B, B such that Rh=Rh(Mascaret) (default)
    @param option_z (int): option for vertical position of rectangular profile
        option_z=1 : zb=min(Zb(y))
        option_z=2 : zb=mean(Zb(y)) (default)
        option_z=3 : zs=Zs
    @param verbose (bool): verbose mode (default: False)
    @param plot_profile (bool): plot profile (default: False)

    @ return x (np.array[nx], float): x abscissa of profiles
    @ return zb_arr (np.array[nx], float): bed elevation
    @ return Bx_arr (np.array[nx], float): width
    @ return zs_arr (np.array[nx], float): free surface elevation
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # assert input
    if zs_in is not None and h_in is not None:
        raise ValueError("set zs_in OR h_in")
    if zs_in is not None and xs_in is None:
        raise ValueError("set xs_in with zs_in")

    # read Mascaret .geo file
    x, nx, ny, profiles = read_mascaret_geo(mascaret_file, verbose=verbose)

    # convert each profile
    zb_arr = np.empty(nx, dtype='d')
    zs_arr = np.empty(nx, dtype='d')
    Bx_arr = np.empty(nx, dtype='d')

    # define interpolation function
    if zs_in is not None and len(xs_in)!=nx:
        f_zs = interp1d(xs_in, zs_in, kind='linear', fill_value="extrapolate")

    for i in range(nx):

        # interpolate zs_in if required
        if zs_in is not None and len(xs_in)!=nx:
            zs_i = f_zs(x[i])
        else:
            zs_i = zs_in[i] if zs_in is not None else None

        # convert profile
        Bx_arr[i], zb_arr[i], _, zs_arr[i] = convert_profile(
            y=profiles[i, 0:ny[i], 0], 
            zb_in=profiles[i, 0:ny[i], 1], 
            zf_in=None, 
            zs_in=zs_i,
            h_in=h_in,
            option=option,
            option_z=option_z,
            verbose=verbose,
            plot_profile=plot_profile)

    # write Arrakis .geo file
    write_geo_file(arrakis_file, x, zb_arr, ZF=None, L=Bx_arr)

    return x, zb_arr, Bx_arr, zs_arr

def convert_geo_courlis_to_arrakis(
        courlis_file,
        arrakis_file="MESH.geo",
        xs_in=None,
        zs_in=None,
        h_in=None,
        option=3,
        option_z=2,
        verbose=False,
        plot_profile=False):
    """ 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Convert Courlis .geoC file into Arrakis .geo file

    Remarks: convertion to rectangular profile requires water elevation data
      - zs_in : array containing water surface elevation for each profile
      or 
      - h_in : same constant water depth for all profiles
      - if None of them is given, zs is set to max(zb)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param courlis_file (str): Courlis geoC file
    @param arrakis_file (str): Arrakis geo file (default: "MESH.geo")
    @param xs_in (np.array[nx], float): x abscissa of zs_in
    @param zs_in (np.array[nx], float): free surface elevation
    @param h_in (float): water depth
    @param option (int): option for convertion into rectangular cross-section
        option=1 : zb=Zb, zs=Zs and B=A/(Zs-Zb)
        option=2 : zb=Zb, B=B(Mascaret) computed where Z_s-Z_b>0
        option=3 : zb=Zb, zs=Zb+A/B, B such that Rh=Rh(Mascaret) (default)
    @param option_z (int): option for vertical position of rectangular profile
        option_z=1 : zb=min(Zb(y))
        option_z=2 : zb=mean(Zb(y)) (default)
        option_z=3 : zs=Zs
    @param verbose (bool): verbose mode (default: False)
    @param plot_profile (bool): plot profile (default: False)

    @return x (np.array[nx], float): x abscissa of profiles
    @return zb_arr (np.array[nx], float): bed elevation
    @return zf_arr (np.array[nx], float): non erodible bed elevation
    @return Bx_arr (np.array[nx], float): width
    @return zs_arr (np.array[nx], float): free surface elevation
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    """
    # assert input
    if zs_in is not None and h_in is not None:
        raise ValueError("set zs_in OR h_in")
    if zs_in is not None and xs_in is None:
        raise ValueError("set xs_in with zs_in")

    # read Courlis .geoC file
    x, nx, ny, profiles = read_courlis_geoC(courlis_file, verbose=verbose)

    # define interpolation function
    if zs_in is not None and len(xs_in)!=nx:
        f_zs = interp1d(xs_in, zs_in, kind='linear', fill_value="extrapolate")

    # convert each profile
    zb_arr = np.empty(nx, dtype='d')
    zf_arr = np.empty(nx, dtype='d')
    zs_arr = np.empty(nx, dtype='d')
    Bx_arr = np.empty(nx, dtype='d')

    for i in range(nx):

        # interpolate zs_in if required
        if zs_in is not None and len(xs_in)!=nx:
            zs_i = f_zs(x[i])
        else:
            zs_i = zs_in[i] if zs_in is not None else None

        # convert profile
        Bx_arr[i], zb_arr[i], zf_arr[i], zs_arr[i] = convert_profile(
            y=profiles[i, 0:ny[i], 0], 
            zb_in=profiles[i, 0:ny[i], 1], 
            zf_in=profiles[i, 0:ny[i], 2],
            zs_in=zs_i,
            h_in=h_in,
            option=option,
            option_z=option_z,
            verbose=verbose,
            plot_profile=plot_profile)

    # write Arrakis .geo file
    write_geo_file(arrakis_file, x, zb_arr, ZF=zf_arr, L=Bx_arr)

    return x, zb_arr, zf_arr, Bx_arr, zs_arr

def refine_geo(geo_file="MESH.geo", new_file=None, dx=0.5, option=1):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Refine mesh
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param geo_file (str): geometry file (default: "MESH.geo")
    @param new_file (str): new geometry file (default: None -> overwrite geo_file)
    @param dx (float): target cell size (default: 0.5)
    @param option (int): 1: keep old mesh nodes
                         2: uniform mesh (don't keep old mesh nodes)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # read geo file
    # ~~~~~~~~~~~~~
    X, ZB, ZF, L = read_geo_file(geo_file)

    # refine mesh
    # ~~~~~~~~~~~
    X_new = refine_mesh(X, dx, option)

    # interpolate data on new mesh
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f_zb = interp1d(X, ZB, kind='linear', fill_value="extrapolate")
    ZB_new = f_zb(X_new)
    if ZF is not None:
        f_zf = interp1d(X, ZF, kind='linear', fill_value="extrapolate")
        ZF_new = f_zf(X_new)
    else:
        ZF_new = None
    if L is not None:
        f_L = interp1d(X, L, kind='linear', fill_value="extrapolate")
        L_new = f_L(X_new)
    else:
        L_new = None

    # write geo file
    # ~~~~~~~~~~~~~~
    if new_file==None:
        new_file = geo_file
    
    write_geo_file(new_file, X_new, ZB_new, ZF=ZF_new, L=L_new)

    return 0

def modify_geo(geo_file="MESH.geo", new_file=None, 
               array_list=["ZB"],
               span_list=None,  
               scaling_list=None,
               translation_list=None):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Modify geo file
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param geo_file (str): geometry file (default: "MESH.geo")
    @param new_file (str): new geometry file (default: None -> overwrite geo_file)
    @param span_list (list of int): smoothing span (default: [3])
    @param array_list (list of str): list of arrays to smooth (default: ["ZB"])
        possible values: "ZB", "ZF", "L"
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    """
    # read geo file
    # ~~~~~~~~~~~~~
    X, ZB, ZF, L = read_geo_file(geo_file)
    ZB_new = ZB
    ZF_new = ZF
    L_new = L

    # smooth data
    # ~~~~~~~~~~~
    if span_list is not None:
        assert len(span_list)==len(array_list)

        if "ZB" in array_list:
            index = array_list.index("ZB")
            span = span_list[index]
            ZB_new = smooth_data(ZB, span)

        if "ZF" in array_list:
            if ZF is None:
                raise ValueError("no ZF in geo file")
            index = array_list.index("ZF")
            span = span_list[index]
            ZF_new = smooth_data(ZF, span)

        if "L" in array_list:
            if L is None:
                raise ValueError("no L in geo file")
            index = array_list.index("L")
            span = span_list[index]
            L_new = smooth_data(L, span)

    # scale data
    # ~~~~~~~~~~
    if scaling_list is not None:
        assert len(scaling_list)==len(array_list)

        if "ZB" in array_list:
            index = array_list.index("ZB")
            scaling = scaling_list[index]
            ZB_new *= scaling

        if "ZF" in array_list:
            if ZF is None:
                raise ValueError("no ZF in geo file")
            index = array_list.index("ZF")
            scaling = scaling_list[index]
            ZF_new *= scaling

        if "L" in array_list:
            index = array_list.index("L")
            scaling = scaling_list[index]
            L_new *= scaling

    # translate data
    # ~~~~~~~~~~~~~~
    if translation_list is not None:
        assert len(translation_list)==len(array_list)

        if "ZB" in array_list:
            index = array_list.index("ZB")
            translation = translation_list[index]
            ZB_new += translation

        if "ZF" in array_list:
            if ZF is None:
                raise ValueError("no ZF in geo file")
            index = array_list.index("ZF")
            translation = translation_list[index]
            ZF_new += translation

        if "L" in array_list:
            index = array_list.index("L")
            translation = translation_list[index]
            L_new += translation

    # write geo file
    # ~~~~~~~~~~~~~~
    if new_file==None:
        new_file = geo_file
    
    write_geo_file(new_file, X, ZB_new, ZF=ZF_new, L=L_new)

    return 0