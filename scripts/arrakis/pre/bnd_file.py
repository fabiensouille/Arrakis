import numpy as np

def build_bnd_file(bnd_file="BC.liq",
                   ntimes=100, t0=0., tf=100.,
                   zs_funct=None,
                   h_funct=None,
                   u_funct=None,
                   q_funct=None,
                   t_funct=None, ntrac=0,
                   qs_funct=None,
                   verbose=False):
    """ 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Build boundary condition file
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param bnd_file (str): boundary file (default: "BC.liq")
    @param ntimes (int): number of times (default: 100)
    @param t0 (int): initial time (default: 0.)
    @param tf (int): final time (default: 100.)
    @param zs_funct (function): free surface function (default: None)
    @param h_funct (function): water depth function (default: None)
    @param u_funct (function): velocity function (default: None)
    @param q_funct (function): flow rate function (default: None)
    @param t_funct (function): tracers function (default: None)
    @param ntrac (int): number of tracers (default: 0)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    # check bnd data type
    # ~~~~~~~~~~~~~~~~~~~
    z_flag = bool(zs_funct is not None)
    h_flag = bool(h_funct is not None)
    u_flag = bool(u_funct is not None)
    q_flag = bool(q_funct is not None)
    t_flag = bool(ntrac > 0)
    qs_flag = bool(qs_funct is not None)

    # bnd type
    if (z_flag or h_flag) and (u_flag or q_flag):
        critical = True
        nfields = 3 + ntrac
    else:
        critical = False
        nfields = 2 + ntrac

    if qs_flag: 
        nfields += 1

    # inititalize labels and units
    field_labels = ["-" for i in range(nfields)]
    field_units  = ["-" for i in range(nfields)]

    # exceptions
    if z_flag and h_flag:
        raise ValueError("set free surface OR water depth")
    if u_funct and q_flag:
        raise ValueError("set velocity OR flow rate")
    if t_flag and ntrac==0:
        raise ValueError("ntrac > 0 required")
    if not t_flag and ntrac>0:
        raise ValueError("ntrac > 0 but no t_funct found")

    # field array [T, Zs or H, U or Q, T(1)...T(ntrac)]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    array = np.empty((nfields, ntimes), dtype='d')

    # Times
    array[0, :] = np.linspace(t0, tf, ntimes)
    field_labels[0] = "T"
    field_units[0]  = "s"

    # Zs or H
    if z_flag:
        field_labels[1] = "ZS"
        field_units[1]  = "m"
        for i in range(ntimes):
            array[1, i] = zs_funct(array[0, i])
    elif h_flag:
        field_labels[1] = "H"
        field_units[1]  = "m"
        for i in range(ntimes):
            array[1, i] = h_funct(array[0, i])

    # U or Q
    if not critical:
        l = 1
    else:
        l = 2
    if u_flag:
        field_labels[l] = "U"
        field_units[l]  = "m/s"
        for i in range(ntimes):
            array[l, i] = u_funct(array[0, i])
    elif q_flag:
        field_labels[l] = "Q"
        field_units[l]  = "m3/s"
        for i in range(ntimes):
            array[l, i] = q_funct(array[0, i])

    # bedload flux
    if qs_flag:
        l += 1
        field_labels[l] = "QS"
        field_units[l]  = "m3/s"
        for i in range(ntimes):
            array[l, i] = qs_funct(array[0, i])

    # tracers
    if t_flag:
        l += 1
        for k in range(ntrac):
            field_labels[l+k] = "T{}".format(k+1)
            field_units[l+k]  = "-"
            for i in range(ntimes):
                array[l+k, i] = t_funct(array[0, i], k)

    if verbose:
        print("fields :", array)
        print("labels :", field_labels)
        print("units  :", field_units)

    # write .liq file
    # ~~~~~~~~~~~~~~~
    f = open(bnd_file, "w")
    f.write("# Liquid boundary file \n")
    # write labels
    for j in range(nfields):
        if j>0: f.write(",")
        f.write("{:}".format(field_labels[j]))
    f.write("\n")
    # write units
    for j in range(nfields):
        if j>0: f.write(",")
        f.write("{:}".format(field_units[j]))
    f.write("\n")
    # write fields
    for i in range(ntimes):
        for j in range(nfields):
            f.write("{:} ".format(array[j, i]))
        f.write("\n")
    f.close()

    return 0
