# plot2d_corr_module.py

def corr2d(corr_x, corr_y, idx):
    
    corr_x = np.array(corr_x[idx])
    corr_y = np.array(corr_y[idx])
    coeff  = np.corrcoef(corr_x.flatten(), corr_y.flatten())
    
    return coeff[0, 1]
    

def scatter_corr(corr_x, corr_y, idx):
    hist2d, xedges, yedges = np.histogram2d(corr_x[idx].flatten(), corr_y[idx].flatten(), bins=50)

    x = 0.5*(xedges[0:-1] + xedges[1:])
    y = 0.5*(yedges[0:-1] + yedges[1:])

    xv, yv = np.meshgrid(x, y)
    
    hist2d = hist2d.T
    count = np.log10(hist2d[hist2d>0]/np.sum(hist2d))
    
    return xv[hist2d>0], yv[hist2d>0], count
    

def plot_corr_Q_Sigma(field, frame, sim_time):
    Toomre_Q             = field.Toomre_Q
    surface_density      = field.surface_density
    surface_density_mean = field.surface_density_mean
    
    corr_x = 1/Toomre_Q
    corr_y = surface_density / surface_density_mean
    idx    = 
    
    coeff  = corr2d(corr_x, corr_y, idx)
    xv, yv, color_scale = scatter_corr(corr_x, corr_y, idx)