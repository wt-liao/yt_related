## def_gamer_cyl_field
## define gamer field in cyl geometry

import numpy as np
import yt.units as u
import yt.utilities.physical_constants as const
from yt.fields.api import ValidateSpatial


### 1.0 fix velocity field
def get_velocity_r(field, data):
    return data['gamer','MomX'] / data['gamer','Dens']

def get_velocity_theta(field, data):
    return data['gamer','MomY'] / data['gamer','Dens']


### 2.0 gradient field
## 2.1 dr field
def get_dvr_dr(field, data):
    fieldData = data['velocity_r']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dr"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[2:,1:-1,1:-1] - fieldData[:-2,1:-1,1:-1]) /
        (2. * data["index", "dr"][1:-1,1:-1,1:-1]))
    
    return new_field 

def get_dvtheta_dr(field, data):
    fieldData = data['velocity_theta']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dr"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[2:,1:-1,1:-1] - fieldData[:-2,1:-1,1:-1]) /
        (2. * data["index", "dr"][1:-1,1:-1,1:-1]))
    
    return new_field   

def get_dAM_dr(field, data):
    fieldData = data['velocity_theta'] * data["index", "r"]
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dr"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[2:,1:-1,1:-1] - fieldData[:-2,1:-1,1:-1]) /
        (2. * data["index", "dr"][1:-1,1:-1,1:-1]))
    
    return new_field   
     
def get_dvz_dr(field, data):
    fieldData = data['velocity_z']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dr"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[2:,1:-1,1:-1] - fieldData[:-2,1:-1,1:-1]) /
        (2. * data["index", "dr"][1:-1,1:-1,1:-1]))
    
    return new_field    


##  2.2 d_theta field
def get_dvr_dtheta(field, data):
    fieldData = data['velocity_r']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dtheta"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[1:-1,2:,1:-1] - fieldData[1:-1,:-2,1:-1]) /
        (2. * data["index", "dtheta"][1:-1,1:-1,1:-1]))
    
    return new_field 

def get_dvtheta_dtheta(field, data):
    fieldData = data['velocity_theta']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dtheta"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[1:-1,2:,1:-1] - fieldData[1:-1,:-2,1:-1]) /
        (2. * data["index", "dtheta"][1:-1,1:-1,1:-1]))
    
    return new_field 

def get_dvz_dtheta(field, data):
    fieldData = data['velocity_z']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dtheta"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[1:-1,2:,1:-1] - fieldData[1:-1,:-2,1:-1]) /
        (2. * data["index", "dtheta"][1:-1,1:-1,1:-1]))
    
    return new_field 

    
##  2.3 dz field
def get_dvr_dz(field, data):
    fieldData = data['velocity_r']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dz"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[1:-1,1:-1,2:] - fieldData[1:-1,1:-1,:-2]) /
        (2. * data["index", "dz"][1:-1,1:-1,1:-1]))
    
    return new_field 

def get_dvtheta_dz(field, data):
    fieldData = data['velocity_theta']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dz"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[1:-1,1:-1,2:] - fieldData[1:-1,1:-1,:-2]) /
        (2. * data["index", "dz"][1:-1,1:-1,1:-1]))
    
    return new_field 

def get_dvz_dz(field, data):
    fieldData = data['velocity_z']
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dz"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[1:-1,1:-1,2:] - fieldData[1:-1,1:-1,:-2]) /
        (2. * data["index", "dz"][1:-1,1:-1,1:-1]))
    
    return new_field 


### 3.0 vorticity field
def get_vorticity_r(field, data):
    # vorticity_r     = dvz_dtheta/radius - dvtheta_dz
    return (data['gas','velocity_z_gradient_theta']/data["index", "r"]) - data['gas','velocity_theta_gradient_z']

def get_vorticity_theta(field, data):
    # vorticity_theta = dvr_dz - dvz_dr
    return data['gas','velocity_r_gradient_z'] - data['gas','velocity_z_gradient_r']

def get_vorticity_z(field, data):
    # vorticity_z     = (dAM_dr - dvr_dtheta)/radius
    return (data['gas','AM_gradient_r'] - data['gas','velocity_r_gradient_theta']) / data["index", "r"]

def get_vorticity(field, data):
    vorticity_r     = data['gas', 'vorticity_r']
    vorticity_theta = data['gas', 'vorticity_theta']
    vorticity_z     = data['gas', 'vorticity_z']
    return np.sqrt(vorticity_r**2.0 + vorticity_theta**2.0 + vorticity_z**2.0)

def get_helicity(feild, data):
    #vorticity_r*vr + vorticity_theta*vtheta + vorticity_z*vz
    return data['vorticity_r']    * data['velocity_r'] + \
           data['vorticity_theta']* data['velocity_theta'] + \
           data['vorticity_z']    * data['velocity_z']

### 4.0 epicycle frequency
def get_kappa_square(field, data):
    fieldData = data['velocity_theta']*data["r"]
    new_field = data.ds.arr(np.zeros(fieldData.shape, dtype=fieldData.dtype),
                            fieldData.units/data["index", "dr"].units)
    new_field[1:-1,1:-1,1:-1] = (
        (fieldData[2:,1:-1,1:-1] - fieldData[:-2,1:-1,1:-1]) /
        (2. * data["index", "dr"][1:-1,1:-1,1:-1]))
    
    return new_field*(2*data["velocity_theta"]/data["r"]**2.0)


def get_kappa(field, data):
    fieldData = data['epicycle_freq_square']
    
    outData = np.sqrt(fieldData)
    outData[fieldData < 0] = 0.0/u.s
    
    return outData