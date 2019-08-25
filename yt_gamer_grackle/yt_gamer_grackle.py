"""
Usage:
fs_yt_gamer_grackle.py -f=<filepath> 
                       -coord=<coordinate> 
                       -op=<optical_depth> 
                       -rate=<three_body_rate>

Note:
works under: gamer cyl coord; 9-species grackle field
"""

from fs_grackle_field import *
from fs_grackle_rad_field import *
from gamer_cyl_field import *
from gamer_energy_field import *
#from gamer_dynamo_field import *


def add_popIII_field(ds, **kwargs):
    optical_depth   = kwargs["optical_depth"]
    three_body_rate = kwargs["three_body_rate"]
    
    add_grackle_field(ds, three_body_rate=three_body_rate)
    add_grackle_rad_field(ds, optical_depth=optical_depth)

def add_cyl_field(ds, **kwargs):
    add_gamer_cyl_field(ds)

def add_energy_field(ds, **kwargs):
    add_gamer_energy_field(ds)
    
def add_dynamo_field(ds, **kwargs):
    add_gamer_dynamo_field(ds)


## main 
def yt_addon_main(ds, **kwargs):
    ## add popIII before cyl; because cyl needs mean molecular weight
    add_popIII_field(ds, **kwargs)
    add_cyl_field(ds, **kwargs)
    add_energy_field(ds, **kwargs)



if __name__ == "__main__":
    
    import yt
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, dest='filepath')
    parser.add_argument('-coord', type=str, dest='coordinate')
    parser.add_argument('-op', type=int, dest='optical_depth')
    parser.add_argument('-rate', type=int, dest='three_body_rate')
    args = parser.parse_args()
    
    filepath        = args.filepath
    coordinate      = args.coordinate
    optical_depth   = args.optical_depth
    three_body_rate = args.three_body_rate
    
    ds = yt.load(filepath)
    yt_addon_main(ds, optical_depth=optical_depth, three_body_rate=three_body_rate)
