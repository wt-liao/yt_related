"""
Usage:
fs_yt_gamer_grackle.py -f=<filepath> -coord=<coordinate> -op=<optical_depth>

Note:
works under: gamer cyl coord; 9-species grackle field
"""

from fs_grackle_field import *
from fs_grackle_rad_field import *


def add_popIII_field(ds, **kwargs):
    optical_depth = kwargs["optical_depth"]
    
    add_grackle_field(ds)
    add_grackle_rad_field(ds, optical_depth = optical_depth)


if __name__ == "__main__":
    
    import yt
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, dest='filepath')
    parser.add_argument('-coord', type=str, dest='coordinate')
    parser.add_argument('-op', type=str, dest='optical_depth')
    args = parser.parse_args()
    
    filepath      = args.filepath
    coordinate    = args.coordinate
    optical_depth = args.optical_depth
    
    ds = yt.load(filepath)
    add_popIII_field(ds, optical_depth = optical_depth)
