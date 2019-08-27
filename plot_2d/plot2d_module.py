## plot2d_module.py

import matplotlib.pyplot as plt
import numpy as np
import yt
import yt.units as u
import yt.utilities.physical_constants as const


def scatter2d(s_x, s_y, idx):
    hist2d, xedges, yedges = np.histogram2d(s_x[idx].flatten(), s_y[idx].flatten(), bins=50)

    x = 0.5*(xedges[0:-1] + xedges[1:])
    y = 0.5*(yedges[0:-1] + yedges[1:])
    xv, yv = np.meshgrid(x, y)

    hist2d = hist2d.T
    count = np.log10(hist2d[hist2d>0]/np.sum(hist2d))
    
    return xv[hist2d>0], yv[hist2d>0], count


def plot_surface_density(field, kw):
    frame, save_dir, sim_time = kw['frame'], kw['save_dir'], kw['sim_time']
    
    radius_2d       = field.radius_2d
    surface_density = field.surface_density
    
    s_x = np.array(radius_2d.in_units('au'))
    s_y = surface_density

    ## 2d profile
    plot_x = np.mean(s_x, axis=1)
    plot_y = np.mean(s_y, axis=1)
    
    ## plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    plt.plot( plot_x, plot_y, c='k', lw=2)
    
    ax.text(0.95, 0.88, r'${\rm Time}=$'+r'$%.1f$'%float(sim_time)+r'${\rm\ yr}$',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, color='k', fontsize=17, fontweight='bold', 
            bbox={'facecolor':'w', 'alpha':0.5, 'pad':3})
    
    plt.xlim(0.5, 60)
    plt.ylim(3e2, 2e4)

    plt.xscale('log')
    plt.yscale('log')

    plt.ylabel(r'$<\Sigma>_\theta$')
    plt.xlabel(r'${\rm r\ (AU)}$')
    
    ## 
    filename = 'SurfaceDensity_2D_'+frame
    plt.tight_layout()
    plt.savefig(save_dir+filename+'.eps', bbox_inches='tight')
    plt.close()
    
    

def plot_ToomreQ(field, idx, kw):
    frame, save_dir, sim_time = kw['frame'], kw['save_dir'], kw['sim_time']
    
    Toomre_Q  = field.Toomre_Q
    radius_2d = field.radius_2d
    
    s_x = np.array(radius_2d.in_units('au'))
    s_x = np.log10(s_x)
    s_y = Toomre_Q
    
    ## setting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ## scatter plot
    xv, yv, color_scale = scatter2d(s_x, s_y, idx)
    plt.scatter(xv, yv, c=color_scale, cmap='viridis', marker='o', s=8 )

    ## 2d profile
    plot_x, plot_y = np.mean(s_x, axis=1), np.mean(s_y, axis=1)
    
    ## annotate
    ax.text(0.95, 0.88, r'${\rm Time}=$'+r'$%.1f$'%float(sim_time)+r'${\rm\ yr}$',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, color='k', fontsize=17, fontweight='bold', 
            bbox={'facecolor':'w', 'alpha':0.5, 'pad':3})

    plt.plot([0, 60], [1.5, 1.5], ls='--', lw=1.5, c='gray')
    plt.plot( plot_x[8:-50], plot_y[8:-50], c='k', ls='-', lw=2)

    plt.ylabel(r'${\rm Toomre\ Q}$')
    plt.xlabel(r'$ \log_{10} r $')

    plt.ylim(0.5, 3)
    plt.xlim(0.25, 1.7)
    
    ## save
    filename = 'ToomreQ_2D_'+frame
    plt.tight_layout()
    plt.savefig(save_dir+filename+'.eps', bbox_inches='tight')
    plt.close()


def plot_kappa(field, idx, kw):
    frame, save_dir, sim_time = kw['frame'], kw['save_dir'], kw['sim_time']
    
    kappa_square_mid  = field.kappa_square_mid
    omega_Kep_mean    = field.omega_Kep_mid
    radius_2d         = field.radius_2d
    
    s_x = np.array(radius_2d.in_units('au'))
    s_x = np.log10(s_x)
    s_y = (kappa_square_mid / omega_Kep_mean**2).in_cgs()
    
    ## setting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ## scatter plot
    xv, yv, color_scale = scatter2d(s_x, s_y, idx)
    plt.scatter(xv, yv, c=color_scale, cmap='viridis', marker='o', s=8 )

    ## 2d profile
    plot_x, plot_y = np.mean(s_x, axis=1), np.mean(s_y, axis=1)
    
    ## annotate
    ax.text(0.42, 0.88, r'${\rm Time}=$'+r'$%.1f$'%float(sim_time)+r'${\rm\ yr}$',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, color='k', fontsize=17, fontweight='bold', 
            bbox={'facecolor':'w', 'alpha':0.5, 'pad':3})

    plt.plot([0, 60], [1, 1], ls='--', lw=1.5, c='gray')
    plt.plot([0, 60], [0, 0], ls=':',  lw=1.5, c='gray')
    
    plt.plot( plot_x[8:-50], plot_y[8:-50], c='k', ls='-', lw=2)

    plt.ylabel(r'$\kappa^2/ {< \Omega_{\rm Kep} >_\theta} ^2$')
    plt.xlabel(r'$ \log_{10} r $')

    plt.ylim(-3.5, 5)
    plt.xlim(0.25, 1.7)
    
    ## save
    filename = 'KappaSquare_'+frame
    plt.tight_layout()
    plt.savefig(save_dir+filename+'.eps', bbox_inches='tight')
    plt.close()
    





