import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial

def add_grackle_field(ds):
    
    mass_unit   = ds.mass_unit.in_cgs()
    length_unit = ds.length_unit.in_cgs()
    time_unit   = ds.time_unit.in_cgs()
    rho_unit    = mass_unit / length_unit**3
    
    N1 = ds.domain_dimensions[0]
    N2 = ds.domain_dimensions[1]
    N3 = ds.domain_dimensions[2]
    
    
    #######################
    #### I. functions #####
    #######################
    # 1.0 fix grackle field
    def get_electron(field, data):
        return data['gamer', 'Electron'] * rho_unit
        
    def get_HI(field, data):
        return data['gamer', 'HI'] * rho_unit
        
    def get_HII(field, data):
        return data['gamer', 'HII'] * rho_unit
    
    def get_HeI(field, data):
        return data['gamer', 'HeI'] * rho_unit
        
    def get_HeII(field, data):
        return data['gamer', 'HeII'] * rho_unit
    
    def get_HeIII(field, data):
        return data['gamer', 'HeIII'] * rho_unit
    
    def get_HM(field, data):
        return data['gamer', 'HM'] * rho_unit
        
    def get_H2I(feild, data):
        return data['gamer', 'H2I'] * rho_unit
    
    def get_H2II(field, data):
        return data['gamer', 'H2II'] * rho_unit
        
    # 2.0 rates
    def get_H2_rate_f(field, data):
        ## H2 formation rate 
        T = data['temperature']/u.K
        
        return 5.5e-29 / T * (u.cm**6/u.s)
        
    def get_H2_rate_d(field, data):
        ## H2 destruction rate
        T = data['temperature']/u.K
        
        return 5.24e-7*T**(-0.485)*np.exp(-52000/T) * (u.cm**3/u.s)
    
    # 3.0 derived chemistry field
    def get_H2_mass_fraction(field, data):
        return data['gas', 'H2I'] / (0.76*data['gas', 'density'])
    
    
    #########################
    #### II. add fields #####
    #########################
    # 1.0 fix grackle field
    ds.add_field(('gas','Electron'), function=get_electron, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','HI'), function=get_HI, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','HII'), function=get_HII, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','HeI'), function=get_HeI, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','HeII'), function=get_HeII, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','HeIII'), function=get_HeIII, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','HM'), function=get_HM, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','H2I'), function=get_H2I, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    ds.add_field(('gas','H2II'), function=get_H2II, \
                 units="g/cm**3", force_override=True, \
                 sampling_type='cell', \
                 validators=[ValidateSpatial(ghost_zones=0)] )
    
    # 2.0 rate
    ds.add_field(('gas', 'H2_formation_rate'), \
                 function=get_H2_rate_f, units="cm**6/s", \
                 sampling_type='cell')
    ds.add_field(('gas', 'H2_destruction_rate'), \
                 function=get_H2_rate_d, units="cm**3/s", \
                 sampling_type='cell')
                 
    # 3.0 derived chemistry field
    ds.add_field(('gas', 'H2_mass_fraction'), \
                 function=get_H2_mass_fraction, units="", take_log=False, \
                 sampling_type='cell')
    
    
    
    