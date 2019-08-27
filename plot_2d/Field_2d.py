## Field_2d.py

import numpy as np
import yt
import yt.units as u
import yt.utilities.physical_constants as const

from yt_gamer_grackle import *

class Field_2d:
    
    def __init__(self, ds, optical_depth=3, three_body_rate=4):
        self.optical_depth   = optical_depth
        self.three_body_rate = three_body_rate
    
        self.N1, self.N2, self.N3 = ds.domain_dimensions
        ds.periodicity = (True, True, True)
    
        yt_addon_main(ds, optical_depth=optical_depth, three_body_rate=three_body_rate)
        
        print('Generate 2D field... ')
        self._init(ds)
        self._mean_field()
        print('Finish 2D field... ')
    
    
    def _init(self, ds):
        self._init1(ds)
        self._init2(ds)
        #self._init3(ds)
    
    
    def _init1(self, ds):
        
        N1, N2, N3 = self.N1, self.N2, self.N3
        ad = ds.covering_grid(0, left_edge=ds.domain_left_edge[:], dims=[N1, N2, N3], num_ghost_zones=1)
        
        # 1..
        dz                    = ad["dz"]
        density               = ad["density"]
        self.surface_density  = np.sum(density*dz, axis=2)
        
        cs                    = ad["sound_speed"]
        kappa                 = ad["epicycle_freq"]
        cs_weighted           = np.sum(cs*density*dz, axis=2)    / self.surface_density
        kappa_weighted        = np.sum(kappa*density*dz, axis=2) / self.surface_density
        
        self.Toomre_Q         = cs_weighted*kappa_weighted/(np.pi*const.G* self.surface_density)
        self.density_mid      = 0.5*( density[:, :, int(0.5*N3)]      + density[:, :, int(0.5*N3)-1] )
        self.kappa_mid        = 0.5*( kappa[:, :, int(0.5*N3)]        + kappa[:, :, int(0.5*N3)-1] )
        
        del dz, density, cs, kappa, cs_weighted, kappa_weighted
        print("\t 1st part of 2D-field done... ")
        
        
    def _init2(self, ds):
        
        N1, N2, N3 = self.N1, self.N2, self.N3
        ad = ds.covering_grid(0, left_edge=ds.domain_left_edge[:], dims=[N1, N2, N3], num_ghost_zones=1)
    
        # 2.
        radius                = ad["r"]
        theta                 = ad["theta"]
        self.theta_2d         = theta[:, :, int(0.5*N3)]
        self.radius_2d        = radius[:, :, int(0.5*N3)]
        del theta
        
        # 3.
        vtheta       = ad["velocity_theta"]
        vtheta_Kep   = ad["velocity_theta_keplerian"]
        omega        = vtheta/radius
        omega_Kep    = vtheta_Kep/radius
        
        self.omega_mid        = 0.5*( omega[:, :, int(0.5*N3)]        + omega[:, :, int(0.5*N3)-1] )
        self.omega_Kep_mid    = 0.5*( omega_Kep[:, :, int(0.5*N3)]    + omega_Kep[:, :, int(0.5*N3)-1] )
        del vtheta, vtheta_Kep, omega, omega_Kep, radius
        
        # 4.
        kappa_square     = ad["epicycle_freq_square"]
        self.kappa_square_mid = 0.5*( kappa_square[:, :, int(0.5*N3)] + kappa_square[:, :, int(0.5*N3)-1] )
        del kappa_square
        
        print("\t 2nd part of 2D-field done... ")
    
    
    def init_rad_field(self, ds):
        
        N1, N2, N3 = self.N1, self.N2, self.N3
        ad = ds.covering_grid(0, left_edge=ds.domain_left_edge[:], dims=[N1, N2, N3], num_ghost_zones=1)
        
        # 3.
        cooling_time     = ad["total_cooling_time"]
        f_esc            = ad["H2_escape_frac_disk"]
        temperature      = ad["temperature"]
        PdV_cooling_time = ad["PdV_cooling_time"]
        
        
        self.cooling_time_mid = 0.5*( cooling_time[:, :, int(0.5*N3)] + cooling_time[:, :, int(0.5*N3)-1] )
        self.f_esc_mid        = 0.5*( f_esc[:, :, int(0.5*N3)]        + f_esc[:, :, int(0.5*N3)-1] )
        self.temperature_mid  = 0.5*( temperature[:, :, int(0.5*N3)]  + temperature[:, :, int(0.5*N3)-1] )
        self.PdV_cooling_time_mid = 0.5*( PdV_cooling_time[:, :, int(0.5*N3)]+ PdV_cooling_time[:, :, int(0.5*N3)-1] )
        
        del cooling_time, f_esc, temperature, PdV_cooling_time
        print("\t 3rd part of 2D-field done... ")
        
        
        
    
    def _mean_field(self):
        surface_density_mean = np.mean(self.surface_density, axis=1)
        density_mean         = np.mean(self.density_mid,     axis=1)
        omega_mean           = np.mean(self.omega_mid,       axis=1)
        omega_Kep_mean       = np.mean(self.omega_Kep_mid,   axis=1)
        
        self.surface_density_mean = np.repeat(surface_density_mean[:, np.newaxis], self.N2, axis=1)
        self.density_mean         = np.repeat(density_mean[:, np.newaxis],         self.N2, axis=1)
        self.omega_mean           = np.repeat(omega_mean[:, np.newaxis],           self.N2, axis=1)      
        self.omega_Kep_mean       = np.repeat(omega_Kep_mean[:, np.newaxis],       self.N2, axis=1)  
        