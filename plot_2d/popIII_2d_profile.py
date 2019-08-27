# 2D profile

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import yt
import glob

import sys
sys.path.insert(0, '/u/sciteam/liao/first_star/python/yt_related/yt_gamer_grackle')

from yt_gamer_grackle import *
from Field_2d import Field_2d
from plot2d_module import *

from pathlib import Path
home = str(Path.home())
#print('Home directory:'+home)


## set matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
    
font = {'weight' : 'bold', 
        'size'   : 18}
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', **font)



## few settings 
def set_save_dir(optical_depth):
    if optical_depth == 3:
        save_dir = home + '/first_star/python/gamer_analysis/popIII/popIII_2dprofile_disk/'
    elif optical_depth == 2:
        save_dir = home + '/first_star/python/gamer_analysis/popIII/popIII_2dprofile_sobolev/'
        
    return save_dir
    
def set_frame_base(optical_depth):
    if optical_depth == 3:
        frame_base = 118
    elif optical_depth == 2:
        frame_base = 57
        
    return frame_base

def get_sim_time(ds):
    current_time = (float(ds.current_time)-0.8)*(159.16319433669406*u.yr)
    
    time = float(current_time.in_units("yr"))
    time = format(time,'.1f')
    
    return time
    

## 
def main(**kwargs):
    base_path       = kwargs["base_path"]
    optical_depth   = kwargs["optical_depth"]
    three_body_rate = kwargs["three_body_rate"]
    
    units_override = {"length_unit": (1.0, "AU"),
                      "time_unit"  : (159.16319433669406, "yr"),
                      "mass_unit"  : (1.0e-6, "msun")}
    
    save_dir = set_save_dir(optical_depth)
        
    filepaths = glob.glob(base_path+'/Data_*')
    
    for path in filepaths:
        print('At file = ' + path)
        ds      = yt.load( path, units_override=units_override )
        field2d = Field_2d(ds, optical_depth=optical_depth, three_body_rate=three_body_rate)
        
        radius_2d, Toomre_Q = field2d.radius_2d, field2d.Toomre_Q
        
        
        idx     = (np.array(radius_2d.in_units('au')) > 2)  & \
                  (np.array(radius_2d.in_units('au')) < 45) & \
                  (Toomre_Q > 0.) & (Toomre_Q < 5.0) 
                  
        frame    = format( int(float(path[-6:]) - set_frame_base(optical_depth)) , '06d')
        sim_time = get_sim_time(ds)
        kw_dict  = {'frame':frame, 'save_dir':save_dir, 'sim_time':sim_time}
        
        ## ploting...
        plot_surface_density(field2d, kw_dict)
        plot_ToomreQ(field2d, idx, kw_dict)
        plot_kappa(field2d, idx, kw_dict)
        
        ## delete obj
        del ds, field2d


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-path', type=str, dest='base_path')
    parser.add_argument('-op', type=int, dest='optical_depth')
    parser.add_argument('-rate', type=int, dest='three_body_rate')
    args = parser.parse_args()
    
    base_path       = args.base_path
    optical_depth   = args.optical_depth
    three_body_rate = args.three_body_rate
        
    main(base_path=base_path, optical_depth=optical_depth, three_body_rate=three_body_rate)
    




