import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial

def add_rad_disk_field(ds):
    mass_unit   = ds.mass_unit.in_cgs()
    length_unit = ds.length_unit.in_cgs()
    time_unit   = ds.time_unit.in_cgs()
    rho_unit    = mass_unit / length_unit**3
    
    #######################
    #### I. functions #####
    #######################
    
    ## 1.0 basic disk rad field:
    def get_disk_op_depth(field, data):
        f_H2 = data['gas', 'H2_mass_fraction']
        Op   = data['gamer', 'DiskOpacity']
        norm = f_H2*float(rho_unit)*float(length_unit)
        
        return Op*norm
    
    def get_H2_disk_escape_frac(field, data):
        tau   = data['gas', 'optical_depth_disk']
        f_esc = (1.0-np.exp(-tau)) / tau
        f_esc[f_esc>1] = 1.0
        
        return f_esc
    
    
    #########################
    #### II. add fields #####
    #########################
    ## double check if simulation is using disk opacity
    #if ( ds.derived_field_list == ('gamer', 'DiskOpacity') ):
    if 1:    
        ds.add_field(('gas', 'optical_depth_disk'), \
                     function=get_disk_op_depth, units="", \
                     sampling_type='cell')
        ds.add_field(('gas', 'H2_escape_frac_disk'), \
                     function=get_H2_disk_escape_frac, units="", \
                     sampling_type='cell')
        
        
        
        
        
        
    
    