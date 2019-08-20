import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial

##
def add_rad_sobolev_field(ds):
    mass_unit   = ds.mass_unit.in_cgs()
    length_unit = ds.length_unit.in_cgs()
    time_unit   = ds.time_unit.in_cgs()
    rho_unit    = mass_unit / length_unit**3
    
    #######################
    #### I. functions #####
    #######################
    
    ## 1.0 basic sobolev field: 
    def get_op_depth_x(field, data):
        n_H2 = data['gas', 'H2I'] / (2.0*const.mass_hydrogen_cgs)
        tau  = data['gamer', 'Opacity_X']
        return np.array(n_H2)*tau*float(length_unit)
    def get_op_depth_y(field, data):
        n_H2 = data['gas', 'H2I'] / (2.0*const.mass_hydrogen_cgs)
        tau  = data['gamer', 'Opacity_Y']
        return np.array(n_H2)*tau*float(length_unit)
    def get_op_depth_z(field, data):
        n_H2 = data['gas', 'H2I'] / (2.0*const.mass_hydrogen_cgs)
        tau  = data['gamer', 'Opacity_Z']
        return np.array(n_H2)*tau*float(length_unit)
    
    ###!! CHECK if the length scale is correct
    def get_sobolev_length_x(field, data):
        tau   = data['gamer', 'Opacity_X']
        alpha = data['gamer', 'Alpha']
        return (tau/alpha*length_unit).in_units('au')
    def get_sobolev_length_y(field, data):
        tau   = data['gamer', 'Opacity_Y']
        alpha = data['gamer', 'Alpha']
        return (tau/alpha*length_unit).in_units('au')
    def get_sobolev_length_z(field, data):
        tau   = data['gamer', 'Opacity_Z']
        alpha = data['gamer', 'Alpha']
        return (tau/alpha*length_unit).in_units('au')
    def get_sobolev_length(field, data):
        lx = data['gas', 'sobolev_length_x']
        ly = data['gas', 'sobolev_length_y']
        lz = data['gas', 'sobolev_length_z']
        return np.sqrt(lx*lx + ly*ly + lz*lz)
        
    def get_H2_sobolev_escape_frac(field, data):
        tau_x = data['gas', 'optical_depth_x']
        tau_y = data['gas', 'optical_depth_y']
        tau_z = data['gas', 'optical_depth_z']
        beta_x = (1-np.exp(-tau_x)) / tau_x
        beta_y = (1-np.exp(-tau_y)) / tau_y
        beta_z = (1-np.exp(-tau_z)) / tau_z
        beta_x[beta_x>1] = 1.0
        beta_y[beta_y>1] = 1.0
        beta_z[beta_z>1] = 1.0
        
        return (beta_x+beta_y+beta_z)/3.0
        
    
    #########################
    #### II. add fields #####
    #########################
    
    ## double check if simulation is using sobolev
    if 1:
        ds.add_field(('gas', 'optical_depth_x'), \
                      function=get_op_depth_x, units="", \
                      sampling_type='cell')
        ds.add_field(('gas', 'optical_depth_y'), \
                      function=get_op_depth_y, units="", \
                      sampling_type='cell')
        ds.add_field(('gas', 'optical_depth_z'), \
                     function=get_op_depth_z, units="", \
                     sampling_type='cell')
                 
        ds.add_field(('gas', 'sobolev_length_x'), \
                     function=get_sobolev_length_x, units="au", \
                     sampling_type='cell')
        ds.add_field(('gas', 'sobolev_length_y'), \
                     function=get_sobolev_length_y, units="au", \
                     sampling_type='cell')
        ds.add_field(('gas', 'sobolev_length_z'), \
                     function=get_sobolev_length_z, units="au", \
                     sampling_type='cell')
        ds.add_field(('gas', 'sobolev_length'), \
                     function=get_sobolev_length, units="au", \
                     sampling_type='cell')
    
        ds.add_field(('gas', 'H2_escape_frac_sobolev'), \
                     function=get_H2_sobolev_escape_frac, units="", \
                     sampling_type='cell')
    
    



