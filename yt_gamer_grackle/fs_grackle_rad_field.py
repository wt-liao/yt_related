import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial

from fs_grackle_rad_sobolev_field import *
from fs_grackle_rad_disk_field import *


def add_grackle_rad_field(ds, optical_depth=1):
    """
    optical_depth: 0:None; 1:fitting function; 2:sobolev; 3:disk
    """
    mass_unit   = ds.mass_unit.in_cgs()
    length_unit = ds.length_unit.in_cgs()
    time_unit   = ds.time_unit.in_cgs()
    rho_unit    = mass_unit / length_unit**3
    
    N1 = ds.domain_dimensions[0]
    N2 = ds.domain_dimensions[1]
    N3 = ds.domain_dimensions[2]
    
    if optical_depth == 2:
        add_rad_sobolev_field(ds)
    elif optical_depth == 3:
        add_rad_disk_field(ds)
    
    
    #######################
    #### I. functions #####
    #######################
    ## 1.2 radiative related
    def get_H2_fitted_escape_frac(field, data):
        n_0 = 8.0e9 /u.cm**3
        n   = ( 0.76*data['gas', 'density']/const.mass_hydrogen_cgs ).in_cgs()
        
        fs = ( (n/n_0).in_cgs() )**(-0.45)        
        fs[fs>1.0] = 1.0
        return fs
    
    
    def get_H2_cooling(field, data):
        rho  = data['gas', 'density'] / (u.g/u.cm**3)
        T    = data['temperature'] / u.K
        f_H2 = data['gas', 'H2_mass_fraction']
        T3   = T*1.0e-3
        
        cool_part  = ( 9.5e-22*T3**(3.76) ) / (1.0+0.12*T3**(2.1)) * np.exp(-(0.13/T3)**3)
        cool_part += 3.0e-24 * np.exp(-0.51/T3) 
        cool_part += 6.7e-19 * np.exp(-5.86/T3)
        cool_part += 1.6e-18 * np.exp(-11.7/T3)
        
        cooling  = 0.76*f_H2/float(const.mass_hydrogen_cgs) * cool_part
        cooling *= rho  ## convert to rate per vol
        
        ## return cooling rate per vol
        unit = u.erg / u.cm**3 / u.s
        return cooling * unit
        
        
    def get_CIE_cooling(field, data):
        rho  = data['gas', 'density'] / (u.g/u.cm**3)
        T    = data['temperature'] / u.K
        f_H2 = data['gas', 'H2_mass_fraction'] 
        
        ## return cooling rate per vol
        unit = u.erg / u.cm**3 / u.s
        return 0.072 * rho**2 * T**4 * f_H2 * 0.76 * unit
    
    
    def get_cooling_rate(field, data):
        cooling_CIE = data['gas', 'CIE_cooling_rate']
        cooling_H2  = data['gas', 'H2_optical_thin_cooling_rate']
        if optical_depth == 2:
            f_esc = data['gas', 'H2_escape_frac_sobolev']
        elif optical_depth == 3:
            f_esc = data['gas', 'H2_escape_frac_disk']

        cooling_tot = f_esc*cooling_H2 + cooling_CIE
        return cooling_tot
    
       
    def get_cooling_time(field, data):
        E_int        = data['gas', 'thermal_energy'] * data['gamer', 'Dens'] # thermal engy per vol
        cooling_rate = data['gas', 'popIII_cooling_rate']
        
        return E_int / cooling_rate
    
    
    #########################
    #### II. add fields #####
    #########################
    
    ds.add_field(('gas', 'H2_optical_thin_cooling_rate'), sampling_type='cell', \
                 function=get_H2_cooling, units="erg/s/cm**3", take_log=True)
    ds.add_field(('gas', 'CIE_cooling_rate'), sampling_type='cell', \
                 function=get_CIE_cooling, units="erg/s/cm**3", take_log=True)
    ds.add_field(('gas', 'popIII_cooling_rate'), sampling_type='cell', \
                 function=get_cooling_rate, units="erg/s/cm**3", take_log=True)
    ds.add_field(('gas', 'popIII_cooling_time'), sampling_type='cell', \
                 function=get_cooling_time, units="s", take_log=True)
    
    
    
    
    
    
    