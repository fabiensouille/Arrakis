import numpy as np
from scipy.interpolate import interp1d

def build_ini_file(geo_file="MESH.geo",
                   ini_file="INIC.ini",
                   zs=0., zs_funct=None, zs_vect=None,  
                   h=1., h_funct=None, h_vect=None, x_vect=None,
                   u=0., u_funct=None, u_vect=None,
                   q=0., q_funct=None, q_vect=None,
                   t=0., t_funct=None, t_vect=None, ntrac=0,
                   friction=0., friction_funct=None, friction_vect=None,
                   time=0., eps=1.e-6, verbose=False):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create initial condition file
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param geo_file (str): geometry file (default: "MESH.geo")
    @param ini_file (str): initial condition file (default: "INI.ini")
    @param x_vect (array): abscissa vector for zs, h, u,
    @param zs (float): constant free surface elevation (default: 0.)
    @param zs_funct (function): free surface function (default: None)
    @param zs_vect (array): free surface vector (default: None)
    @param h (float): constant water depth (default: 0.)
    @param h_funct (function): water depth function (default: None)
    @param h_vect (array): water depth vector (default: None)
    @param u (float): constant velocity (default: 0.)
    @param u_funct (function): velocity function (default: None)
    @param u_vect (array): velocity vector (default: None)
    @param q (float): constant flow rate (default: 0.)
    @param q_funct (function): flow rate function (default: None)
    @param q_vect (array): flow rate vector (default: None)
    @param t (float): constant tracers (default: 0.)
    @param t_funct (function): tracers function (default: None)
    @param t_vect (2D array): tracers vector (default: None)
    @param friction (float): constant friction coefficient (default: 0.)
    @param friction_funct (function): friction coefficient function (default: None)
    @param friction_vect (array): friction coefficient vector (default: None)
    @param ntrac (int): number of tracers (default: 0)
    @param time (function): initial time (default: 0.)
    @param eps (float): epsilon for u=q/h (default: 1.e-6)
    @param verbose (bool): verbose mode (default: False)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # check wether water depth or free surface is given
    if zs!=0. or zs_funct is not None or zs_vect is not None: 
        zs_flag = True
    else:
        zs_flag = False

    if zs!=0. and h!=0:
        raise ValueError("Give either zs or h")
    if zs_funct is not None and h_funct is not None:
        raise ValueError("Give either zs_funct or h_funct")
    if zs_vect is not None and h_vect is not None:
        raise ValueError("Give either zs_vect or h_vect")
        
    # check wether velocity or flowrate is given
    if  u!=0. or u_funct is not None or u_vect is not None: 
        vel_flag = True
    else:
        vel_flag = False

    if u!=0. and q!=0.:
        raise ValueError("Give either u or q")
    if u_funct is not None and q_funct is not None:
        raise ValueError("Give either u_funct or q_funct")
    if u_vect is not None and q_vect is not None:
        raise ValueError("Give either u_vect or q_vect")
    
    # read variable in file
    f = open(geo_file, "r")
    f.readline()
    var = f.readline().strip().split(",")
    f.close()
    if 'X' not in var:
        raise ValueError("X not in mesh file")

    # read geometry file abscissa
    data = np.loadtxt(geo_file, skiprows=2)
    x  = data[:, 0]
    nx = len(x)

    # check h, u, q, t vectors and create interpolating functions if needed
    if zs_vect is not None and len(zs_vect)!=nx:  
        if x_vect is None:
            raise ValueError("x_vect required if zs_vect given and len(zs_vect)!=nx")
        fzs = interp1d(x_vect, zs_vect, kind='linear', fill_value="extrapolate")
    else:
        fzs = None

    if h_vect is not None and len(h_vect)!=nx:
        if x_vect is not None:
            raise ValueError("x_vect required if h_vect given and len(h_vect)!=nx")
        fh = interp1d(x_vect, h_vect, kind='linear', fill_value="extrapolate")
    else:
        fh = None
        
    if u_vect is not None and len(u_vect)!=nx:  
        if x_vect is None:
            raise ValueError("x_vect required if u_vect given and len(u_vect)!=nx")
        fu = interp1d(x_vect, u_vect, kind='linear', fill_value="extrapolate")  
    else:
        fu = None
        
    if q_vect is not None and len(q_vect)!=nx:  
        if x_vect is None:
            raise ValueError("x_vect required if q_vect given and len(q_vect)!=nx")
        fq = interp1d(x_vect, q_vect, kind='linear', fill_value="extrapolate")
    else:
        fq = None
        
    if t_vect is not None and np.shape(t_vect) != (ntrac, nx):  
        if x_vect is None:
            raise ValueError("x_vect required if t_vect given and shape(t_vect)!=(ntrac,nx)")
        ft = [interp1d(x_vect, t_vect[k, :], 
                       kind='linear', 
                       fill_value="extrapolate") for k in range(ntrac)]
    else:
        ft = None

    if friction_vect is not None and len(friction_vect)!=nx:  
        if x_vect is None:
            raise ValueError("x_vect required if friction_vect given and len(friction_vect)!=nx")
        ff = interp1d(x_vect, friction_vect, kind='linear', fill_value="extrapolate")
    else:
        ff = None

    # create initial conditions
    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    friction_flag = bool(friction!=0. or friction_funct or friction_vect)
    if friction_flag:
        nfields = 3 + ntrac
    else:
        nfields = 2 + ntrac

    array = np.empty((nfields, nx), dtype='d')
    field_labels = ['' for i in range(nfields)]

    for i in range(nx):
        l = 0

        # init water depth
        if zs_flag:
            field_labels[l] = "ZS"
            if zs_funct is not None:
                array[l, i] = zs_funct(x[i])
            elif zs_vect is not None:
                if fzs is not None:
                    array[l, i] = fzs(x[i])
                else:
                    array[l, i] = zs_vect[i]
            else:
                array[l, i] = zs
        else:
            field_labels[l] = "H"
            if h_funct is not None:
                array[l, i] = h_funct(x[i])
            elif h_vect is not None:
                if fh is not None:
                    array[l, i] = fh(x[i])
                else:
                    array[l, i] = h_vect[i]
            else:
                array[l, i] = h

        l += 1

        # flow rate or velocity
        if vel_flag:
            # init velocity
            field_labels[l] = "U"
            if u_funct is not None:
                array[l, i] = u_funct(x[i])
            elif u_vect is not None:
                if fu is not None:
                    array[l, i] = fu(x[i])
                else:
                    array[l, i] = u_vect[i]
            else:
                array[l, i] = u
        else:
            # init flow rate
            field_labels[l] = "Q"
            if q_funct is not None:
                array[l, i] = q_funct(x[i])
            elif q_vect is not None:
                if fq is not None:
                    array[l, i] = fq(x[i])
                else:
                    array[l, i] = q_vect[i]
            else:
                array[l, i] = q
        l += 1

        # init water depth
        if friction_flag:
            field_labels[l] = "FC"
            if friction_funct is not None:
                array[l, i] = friction_funct(x[i])
            elif friction_vect is not None:
                if ff is not None:
                    array[l, i] = ff(x[i])
                else:
                    array[l, i] = friction_vect[i]
            else:
                array[l, i] = friction
            l += 1

        # tracers
        if ntrac > 0:
            for k in range(ntrac):
                field_labels[l+k] = "T{}".format(k+1)
                if t_funct is not None:
                    array[l+k, i] = t_funct(x[i], k)
                elif t_vect is not None:
                    if ft is not None:
                        array[l+k, i] = ft[k](x[i])
                    else:
                        array[l+k, i] = t_vect[k, i]
                else:
                    array[l+k, i] = t

    if verbose:
        print("fields :", array)
        print("labels :", field_labels)

    # write ini file
    # ~~~~~~~~~~~~~~
    f = open(ini_file, "w")
    f.write("{}, {} \n".format(time, nx))
    # write labels
    f.write("X,")
    for j in range(nfields):
        f.write("{:}".format(field_labels[j]))
        if j<nfields-1:
            f.write(",")
        else:
            f.write("\n")
    # write fields
    for i in range(nx):
        f.write("{:.8f}, ".format(x[i]))
        for j in range(nfields):
            f.write("{:}".format(array[j, i]))
            if j<nfields-1:
                f.write(", ")
            else:
                f.write("\n")

    return 0