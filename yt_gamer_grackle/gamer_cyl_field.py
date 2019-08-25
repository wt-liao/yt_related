import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial

from def_gamer_cyl_field import *

def add_gamer_cyl_field(ds):
    ## force periodicity to manipulate boundary
    ds.periodicity = (True, True, True)
    
    ### 1.0 add velocity field
    ds.add_field(('gas','velocity_r'), function=get_velocity_r, units="cm/s", force_override=False,
                 sampling_type='cell' )
    ds.add_field(('gas','velocity_theta'), function=get_velocity_theta, units="cm/s", force_override=False,
                 sampling_type='cell' )
    
    ### 2.0 velocity gradient field
    ds.add_field(('gas','velocity_r_gradient_r'), function=get_dvr_dr, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','r_velocity_r_gradient_r'), function=get_drvr_dr, units="cm/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_theta_gradient_r'), function=get_dvtheta_dr, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','AM_gradient_r'), function=get_dAM_dr, units="cm/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_z_gradient_r'), function=get_dvz_dr, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_r_gradient_theta'), function=get_dvr_dtheta, units="cm/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_theta_gradient_theta'), function=get_dvtheta_dtheta, units="cm/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_z_gradient_theta'), function=get_dvz_dtheta, units="cm/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_r_gradient_z'), function=get_dvr_dz, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_theta_gradient_z'), function=get_dvtheta_dz, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_z_gradient_z'), function=get_dvz_dz, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ### 3.0 vorticity field
    ds.add_field(('gas','vorticity_r'), function=get_vorticity_r, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','vorticity_theta'), function=get_vorticity_theta, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','vorticity_z'), function=get_vorticity_z, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','vorticity'), function=get_vorticity, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','velocity_divergence'), function=get_velocity_divergence, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    ds.add_field(('gas','helicity'), function=get_helicity, units="cm/s**2", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)] )
    
    
    ### 4.0 epicycle field
    ds.add_field(('gas','epicycle_freq_square'), function=get_kappa_square, units="1/s**2", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)])
    
    ds.add_field(('gas','epicycle_freq'), function=get_kappa, units="1/s", force_override=True,
                 sampling_type='cell',
                 validators=[ValidateSpatial(ghost_zones=1)])
    
               
    ### 5.0 keplerian velocity
    ds.add_field(('gas','velocity_theta_keplerian'), function=get_keplerian_velocity, units="cm/s", force_override=True,
                 sampling_type='cell', take_log=True,
                 validators=[ValidateSpatial(ghost_zones=1)])
    
    
    ### 6.0 energy related field
    ds.add_field(('gas','PdV_cooling'), function=get_PdV_cooling, units="erg/s/cm**3",
                 sampling_type='cell', take_log=True, force_override=True,
                 validators=[ValidateSpatial(ghost_zones=1)])
    
    ds.add_field(('gas','PdV_cooling_time'), function=get_PdV_cooling_time, units="s",
                 sampling_type='cell', take_log=True, force_override=True,
                 validators=[ValidateSpatial(ghost_zones=1)])
    
    