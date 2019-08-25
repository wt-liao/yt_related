import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial

def add_gamer_energy_field(ds):
    
    def get_total_cooling_rate(field, data):
        PdV_cooling     = data['gas', 'PdV_cooling']
        H2_diss_cooling = data['gas', 'H2_dissociation_cooling']
        popIII_cooling  = data['gas', 'popIII_cooling_rate']
        
        return PdV_cooling + H2_diss_cooling + popIII_cooling
        
    def get_total_cooling_time(field, data):
        E_int   = data['gas', 'thermal_energy'] * data['gas', 'density'] # thermal engy per vol
        cooling = data['gas', 'total_cooling_rate']
        
        return E_int / cooling
    
        
    
    ### energy related field
    ds.add_field(('gas','total_cooling_rate'), function=get_total_cooling_rate, 
                 units="erg/s/cm**3", force_override=True,
                 sampling_type='cell', take_log=True,
                 validators=[ValidateSpatial(ghost_zones=1)])
    
    ds.add_field(('gas','total_cooling_time'), function=get_total_cooling_time, 
                 units="s", force_override=True,
                 sampling_type='cell', take_log=True,
                 validators=[ValidateSpatial(ghost_zones=1)])