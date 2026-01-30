#!/usr/bin/python3
'''
Define paths to use in all files when run from laptop, harddrive, or borah.
'''
import os

path_home = None
path_processed = None
path_mseed = None
path_coords = None
path_output = None
path_figures = None
path_gis_dem = None
path_gis_unit = None

def set_paths(location):
    global path_home, path_processed, path_mseed, path_coords, path_output, path_figures, path_gis_dem, path_gis_unit

    if location == 'laptop':
        path_home = os.path.join("/", "home", "mad", "Documents", "research", "wenas-prescribed-burn")
    elif location == 'borah':
        path_home = os.path.join("/", "bsuhome", "madelinehunt", "wenas-prescribed-burn")
    elif location == 'harddrive':
        path_home = os.path.join("/", "media", "mad", "LaCie 2 LT", "research", "wenas-prescribed-burn")
    else:
        raise ValueError(f"Unknown location: {location}")

    path_processed = os.path.join(path_home, "data", "processed")
    path_mseed = os.path.join(path_home, "data", "mseed")
    path_coords = os.path.join(path_home, "data", "gps", "wenas_burn_coords.csv")
    #path_output = os.path.join(path_home, "data", "output")
    path_figures = os.path.join(path_home, "figures")
    path_gis_dem = os.path.join(path_home, "gis", "dem", "USGS_13_n47w121_20250813.tif")
    path_gis_unit = os.path.join(path_home, "gis", "burn_perimeter", "wenas_perimeter.shp")


