# -*- coding: utf-8 -*-
"""
Post-processing functions 

Author: I60179
"""
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc

def set_rcparams(use_latex=True):
    """ Set plot style parameters """
    rc('figure', autolayout=True)
    if use_latex:
        rc('text.latex', preamble=r'\usepackage{lmodern}')
        rc('text', usetex=True)
    rc('font', size=15)
    rc('legend', fontsize=11)
    rc('legend', framealpha=0.25)
    rc('axes', xmargin=0.)
    rc('axes', ymargin=0.)
    rc('axes', labelsize=15)
    rc('axes', labelpad=4.0)
    rc('grid', color='b0b0b0')
    rc('grid', linestyle=':')
    rc('grid', linewidth=.8)
    rc('grid', alpha=1.)

def get_default_color_palette():
    colors1 = ['#002d74', '#e85113', '#1fa12e', '#c9d200', '#f49e00', '#006ab3', '#3381ff', '#f2855a', '#54de64']
    colors2 = ['#3381ff', '#f2855a', '#54de64', '#f5ff33', '#ffb833', '#33adff', '#99c0ff', '#f7b9a1', '#a9efb1']
    colors3 = ['#99c0ff', '#f7b9a1', '#a9efb1', '#faff99', '#ffdb99', '#99d6ff', '#f7b9a1', '#a9efb1', '#faff99']
    return colors1, colors2, colors3

def get_default_markers():
    markers = ['o', '^', 's', 'D', 'v', '<', '>', 'd', 'H']
    return markers

def plot_res(time, res, var, figname, 
             ref=None,
             show=False,
             plot_sediments=False,
             fancy_plot=True,
             annotate_time=False,
             ylim_h=None, 
             ylim_u=None, 
             dpi=300):

    # extract values
    x = None  
    h = None
    u = None
    zb= None
    zf= None
    qs= None
    taub= None
    taus= None

    nvar = len(var)
    
    for k, var in enumerate(var):
        if var=='X':
            x = res[:,k]
            nx = len(x)
        if var=='H':
            h = res[:,k]
        if var=='U':
            u = res[:,k]
        if var=='ZB':
            zb= res[:,k]
        if var=='ZF' and plot_sediments:
            zf =res[:,k]
        if var=='QS' and plot_sediments:
            qs =res[:,k]
        if var=='TB' and plot_sediments:
            taub =res[:,k]
        if var=='TS' and plot_sediments:
            taus =res[:,k]

        if ref is not None:
            if var=='X':
                x_ref = ref[k]
            if var=='H':
                h_ref = ref[k]
            if var=='U':
                u_ref = ref[k]
            if var=='ZB':
                zb_ref= ref[k]

    # check extracted values
    if x is None:
        raise ValueError('Missing X in result file')
    if h is None:
        raise ValueError('Missing H in result file')
    if u is None:
        raise ValueError('Missing U in result file')
    if zb is None:
        print('WARNING: Missing ZB in result file')
        zb = np.zeros(nx)
    if zf is None:
        if plot_sediments:
            print('WARNING: Missing ZF in result file')
        if zb is not None:
            zf = zb
        else:
            zf = np.zeros(nx)

    fr = [abs(u[i])/np.sqrt(9.80665*max(h[i], 1e-8)) for i in range(len(h))]

    # plot layout
    if plot_sediments:
        ngraph = 3
        figsize = (8, 8)
    else:
        ngraph = 2
        figsize = (8, 6)

    set_rcparams()
    colors, _, _ = get_default_color_palette()
    fig, axes = plt.subplots(ngraph, 1, figsize=figsize)
    ax1 = axes[0]
    ax2 = axes[1]

    # plot h
    if ylim_h is not None:
        zmin = min(min(zb), min(zf), ylim_h[0])
    else:
        zmin = min(min(zb), min(zf))

    if plot_sediments:
        ax1.plot(x, zb, label='$z_b$', color='peru', lw=0.5)
        ax1.plot(x, zf, label='$z_f$', color='saddlebrown', lw=0.5)
        if fancy_plot:
            ax1.fill_between(x, zb, zf, color='peru', alpha=0.45)
            ax1.fill_between(x, zf, zmin, color='saddlebrown', alpha=0.65)
    else:
        ax1.plot(x, zb, label='$z_b$', color='k', lw=0.5)
        if fancy_plot:
            ax1.fill_between(x, zb, zmin, color='grey', alpha=0.5)

    if ref is not None:
        ax1.plot(x_ref, zb_ref+h_ref, label='$z_b+h$ (ref)', color=colors[0], ls=':', lw=1.)

    ax1.plot(x, zb+h, label='$z_b+h$', color=colors[0], marker='o', markersize=1, lw=1.)
    if fancy_plot:
        ax1.fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)

    ax1.legend()
    ax1.set_ylabel("$h$ (m)")
    ax1.set_xlabel("$x$ (m)")
    if ylim_h is not None:
        ax1.set_ylim([ylim_h[0], ylim_h[1]])
    else:
        if plot_sediments:
            ax1.set_ylim([min(zb-0.5*h), max(zb+h)])
        else:
            ax1.set_ylim([min(zb), max(zb+h)])
    ax1.grid()

    # annotate time
    if annotate_time:
        xlim_h = [np.min(x), np.max(x)]
        if ylim_h is None:
            ylim_h = [np.min(zb), np.max(zb+h)]
        
        x_text = xlim_h[0]+0.75*(xlim_h[1]-xlim_h[0])
        y_text = ylim_h[0]+0.85*(ylim_h[1]-ylim_h[0])
        ax1.text(x_text, y_text, '$Time={:10.2}$s'.format(time), fontdict=None)

    # plot u
    if ref is not None:
        ax2.plot(x_ref, h_ref*u_ref, label='$q$ (ref)', color=colors[1], ls=':')
        ax2.plot(x_ref, u_ref, label='$u$ (ref)', color=colors[2], ls=':')

    ax2.plot(x, h*u, label='$q$', color=colors[1], marker='o', markersize=1, lw=1.)
    ax2.plot(x, u, label='$u$', color=colors[2], marker='s', markersize=1, lw=1.)
    ax2.legend(loc=2)
    ax2.set_ylabel("$q$ (m$^2$/s), $u$ (m/s)")
    ax2.set_xlabel("$x$ (m)")
    if ylim_u is not None:
        ax1.set_ylim=(ylim_u[0], ylim_u[1])
    ax2.grid()

    ax2t = ax2.twinx()
    ax2t.set_ylabel("$Fr$ (-)")
    ax2t.plot(x, fr, label='$Fr$', color=colors[3], marker='^', markersize=1, lw=1.)
    ax2t.legend()

    # plot sediment parameters
    if plot_sediments:
        ax3 = axes[2]
        if qs is not None:
            ax3.plot(x, qs, label='$q_s$', color=colors[1], marker='o', markersize=1, lw=1.)
            ax3.set_ylabel("$q_s$ (m$^2$/s)")
            ax3.legend(loc=2)

        if taus is not None:
            ax3t = ax3.twinx()
            ax3t.plot(x, taus, label='$\\tau^*$', color=colors[3], marker='x', markersize=1, lw=1.)
            ax3t.set_ylabel("$\\tau^*$ (-)")
            ax3t.legend(loc=1)
        else:
            if taub is not None:
                ax3t = ax3.twinx()
                ax3t.plot(x, taub, label='$\\tau_b$', color=colors[3], marker='x', markersize=1, lw=1.)
                ax3t.set_ylabel("$\\tau_b$ (Pa)")
                ax3t.legend(loc=1)

        ax3.set_xlabel("$x$ (m)")
        ax3.grid()

    # save figure
    plt.savefig(figname, dpi=dpi)
    if show:
        plt.show()

    # Close figure:
    plt.close('all')


def plot_trac(reshydro, restrac, figname, ref=None,
             figsize=(12, 6), show=False, 
             ylim_h=None, ylim_u=None, dpi=300):
             
    # extract values
    x = reshydro[:,0]
    zb= reshydro[:,1]
    h = reshydro[:,2]
    t = restrac[:,1]
    if ref is not None:
        x_ref = ref[0]
        T_ref = ref[1]
        
    # plot layout
    set_rcparams()
    colors, _, _ = get_default_color_palette()
    fig, axes = plt.subplots(2, 1, figsize=figsize)
    ax1 = axes[0]
    ax2 = axes[1]
    
    # plot h
    ax1.plot(x, zb,label='z', color='k')
    ax1.plot(x, zb+h, label='z+h', color=colors[0], marker='+')
    ax1.legend(loc=2)
    ax1.set_ylabel("h (m)")
    ax1.set_xlabel("x (m)")
    if ylim_h is not None:
      ax1.set_ylim=(ylim_h[0], ylim_h[1])
      
    # plot t
    if ref is not None:
        ax2.plot(x_ref, T_ref, label='T (ref)', color=colors[1], ls=':')
    ax2.plot(x, t, label='T', color='b', marker='+')
    ax2.legend(loc=2)
    ax2.set_ylabel("trac")
    ax2.set_xlabel("x (m)")
    if ylim_h is not None:
      ax2.set_ylim=(ylim_h[0], ylim_h[1])
      
    # save figure
    plt.savefig(figname, dpi=dpi)
    if show:
        plt.show()
    plt.close()

def plot_convergence(dx, error, label='',
                     xlim=None, ylim=None, figsize=(6, 4),
                     xlabel="$\Delta x$",
                     ylabel="Error",
                     plot_ref_order=True,
                     figname='fig.png',
                     dpi=300):

    set_rcparams()
    colors, _, _ = get_default_color_palette()
    markers = get_default_markers()
    
    # plot errors:
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    # plot
    if isinstance(error, list):
        nx = len(dx[0])
        dx_ref = dx[0]
        err_ref1 = error[0][0]
        err_ref2 = error[0][1]
        for i, err in enumerate(error):
            ax.plot(dx[i], err, label=label[i], color=colors[i], marker=markers[i], markersize=3)
    else:
        nx = len(dx)
        err_ref1 = error[0]
        err_ref2 = error[1]
        ax.plot(dx, error, label=label, color=colors[0], marker=markers[0], markersize=3)

    # plot order 1 and 2 ref curves
    if plot_ref_order:
        ord1 = [err_ref1*dx_ref[i]/dx_ref[0] for i in range(nx)]
        ax.plot(dx_ref, ord1, label='order 1', color='r', ls='--', lw=0.5, marker='')
        ord2 = [err_ref2*(dx_ref[i]/dx_ref[0])**2 for i in range(nx)]
        ax.plot(dx_ref, ord2, label='order 2', color='b', ls='--', lw=0.5, marker='')

    # Axe labels
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
        
    # fancy grid
    plt.grid(which='major', color='grey', linestyle='--')
    plt.grid(which='minor', color='grey', linestyle=':')
    
    # Setting x,y scale
    plt.xscale('log')
    plt.yscale('log')
    
    # Legend
    plt.legend()
    # Save figure
    plt.savefig(figname, dpi=dpi)
    plt.show()
    
    # Close figure:
    plt.close('all')


def plotbar(
        data, fig_size=None,
        fig_name='', fig_title=None,
        x_labels='', y_label='', ylim=None,
        legend_labels=None, split_bars=True,
        bar_width=.75,
        y_scale='linear',
        annotate=False,
        annotate_format='e',
        annotate_threshold=None,
        **kwargs):
    """
    Plot data with bars
    """
    set_rcparams()
    colors, colors2, colors3 = get_default_color_palette()
    
    # plot initialization
    fig, ax = plt.subplots(1, 1, figsize=fig_size)

    # number of pos and bars per pos
    if isinstance(data[0], list):
        npos = len(data[0])
        nbar = len(data)
    else:
        npos = len(data)
        nbar = 1

    # plot
    pos = np.arange(npos, dtype='float64')
    barps = []
    for n in range(nbar):
        # split bars on position
        if split_bars:
            d_x = bar_width/nbar
            if nbar % 2 == 0:
                posn = pos + float(n - nbar//2)*d_x + d_x/2
            else:
                posn = pos + float(n - nbar//2)*d_x
        # plot all bars at same position
        else:
            posn = pos
            d_x = bar_width

        barp = plt.bar(posn, data[n], d_x, **kwargs)
        barps.append(barp[0])

        # annotate bars with corresponding values
        if annotate:
            for x, y in zip(posn, data[n]):
                if annotate_format == 'e':
                    label = "{:.2e}".format(y)
                if annotate_format == 'f':
                    label = "{:.2f}".format(y)
                if annotate_format == 'f1':
                    label = "{:.1f}".format(y)
                if annotate_format == 'f0':
                    label = "{:.0f}".format(y)
                if annotate_threshold is not None:
                    if y <= annotate_threshold:
                        plt.annotate(
                            label, (x, y), textcoords="offset points",
                            fontsize=10, xytext=(0, 5), ha='center')
                else:
                    plt.annotate(
                        label, (x, y), textcoords="offset points",
                        fontsize=10, xytext=(0, 5), ha='center')

    # Setting y scale
    plt.yscale(y_scale)
    if y_scale == 'log':
        plt.grid(which='major', color='grey', linestyle='--')
        plt.grid(which='minor', color='grey', linestyle=':')

    # plot options
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])

    # labels
    if y_label != '':
        plt.ylabel(y_label)
    if x_labels != '':
        plt.xticks(pos, x_labels)
    if legend_labels is not None:
        plt.legend((barps), legend_labels)

    if fig_title is not None:
        ax.set_title(fig_title)

    plt.xticks(rotation = 90)
    plt.yticks(rotation = 90)

    # save figure:
    if fig_name != '':
        fig.savefig(fig_name)
    else:
        plt.show()

    # Close figure:
    plt.close('all')

