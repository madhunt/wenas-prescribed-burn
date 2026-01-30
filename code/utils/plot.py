#!/usr/bin/python3
'''Utility functions for plotting.'''

import obspy, os, rasterio, pyproj
import numpy as np
#import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import obspy.signal.array_analysis
import rasterio.warp


def plot_array_response(coords, flim, ax):
    '''
    Calculate and plot the array response on provided axes (with wavenumber in rad/km). 
    Uses sound speed of 343 m/s.
    INPUTS
        coords  : pandas df : Dataframe containing stations for one array/network. 
            Should contain columns for "Latitude" and "Longitude".
        flim    : float     : Frequency limit for bounds of response in Hz.
        ax      : pyplot ax : Axes to plot array response. Plot will show semblance 
            (colormap) across wavenumber in rad/km (x/y axes).
    RETURNS
        res     : pyplot    : Pyplot imshow of array response.
    '''
    # get lat/long coordinates of each station in correct format
    resp_coords = coords[["Longitude", "Latitude"]].to_numpy()
    resp_coords = np.insert(resp_coords, 2, 0, axis=1)

    # calculate wavenumber limit (in rad/km)
    klim = 2 * np.pi * flim / 0.343   
    kstep = 1
    # calculate array response
    array_response = obspy.signal.array_analysis.array_transff_wavenumber(resp_coords,
                                                                            klim=klim,
                                                                            kstep=kstep,
                                                                            coordsys='lonlat')
    # plot array response
    res = ax.imshow(array_response.transpose(), 
                        extent=[-klim, klim, -klim, klim],
                        origin='lower',
                        vmin=0, vmax=1,
                        cmap='plasma')
    return res






def read_dem(path_dem, bounds_region):
    '''
    Reads in DEM from .tif file.
    INPUTS
        path_dem        : str           : Path to DEM .tif file that includes region of interest.
        bounds_region   : list of float : UTM boundary of region of interest within DEM. In UTM 
            coordinates as [min_easting, max_easting, min_northing, max_northing].
    RETURNS
        easting         : np array [n_north, n_east]: Meshgrid of easting coordinates in UTM where elev is sampled.
        northing        : np array [n_north, n_east]: Meshgrid of northing coordinates in UTM where elev is sampled.
        elev            : np array [n_north, n_east]: DEM elevation cropped to region of interest and reprojected in 
            UTM coordinate reference system.
        target_crs      : str           : Final Coordinate Reference System (CRS) of the region in UTM.
    '''
    with rasterio.open(path_dem) as src:
        # get coordinate reference system (CRS) in UTM coords
        bounds_dem = src.bounds
        center_lon = (bounds_dem.left + bounds_dem.right) / 2
        center_lat = (bounds_dem.bottom + bounds_dem.top) / 2
        utm_zone = int((center_lon+180)/6)+1
        hemi = 'north' if center_lat >=0 else 'south'
        target_crs = f"EPSG:{32600+utm_zone if hemi=='north' else 32700+utm_zone}"

        # transform bounds of region of interest from UTM to original CRS of DEM
        transformer = pyproj.Transformer.from_crs(target_crs, src.crs, always_xy=True)
        min_lon, max_lat = transformer.transform(bounds_region[0], bounds_region[3])
        max_lon, min_lat = transformer.transform(bounds_region[1], bounds_region[2])

        # create a window around region of interest in original CRS
        window = rasterio.windows.from_bounds(left=min_lon, bottom=min_lat, 
                                              right=max_lon, top=max_lat, 
                                              transform=src.transform)
        # crop DEM to region of interest and read in
        elev_crop = src.read(1, window=window)
        transform_crop = src.window_transform(window)

        # get bounds of cropped region and calc new reprojection transform
        bounds_crop = rasterio.windows.bounds(window, src.transform)
        transform_utm, width, height = rasterio.warp.calculate_default_transform(src.crs, target_crs, 
                                                                                 elev_crop.shape[1], elev_crop.shape[0], 
                                                                                 *bounds_crop)
        # reproject cropped DEM to UTM CRS and get elevation
        elev = np.empty((height, width), dtype=src.dtypes[0])
        rasterio.warp.reproject(source=elev_crop, 
                                destination=elev, 
                                src_transform=transform_crop, 
                                src_crs=src.crs, 
                                dst_transform=transform_utm, 
                                dst_crs=target_crs, 
                                resampling=rasterio.warp.Resampling.bilinear)

        # get meshgrid of easting and northing in UTM
        x_coords = np.arange(width) * transform_utm.a + transform_utm.c
        y_coords = np.arange(height) * transform_utm.e + transform_utm.f
        easting, northing = np.meshgrid(x_coords, y_coords)
        return easting, northing, elev, target_crs


def plot_contours_from_dem(ax, easting_dem, northing_dem, elev_dem, contour_int=20):
    '''
    Plots contour intervals from DEM on specified figure axis. 
    INPUTS
        ax              : pyplot axis               : Matplotlib axes to plot contour intervals on.
        easting_dem     : np array [n_north, n_east]: Meshgrid of easting coordinates of DEM. Each column 
            contains n_north elements of the same easting coordinate.
        northing_dem    : np array [n_north, n_east]: Meshgrid of northing coordinates of DEM.
        elev_dem        : np array [n_north, n_east]: Elevation of DEM measured at each coordinate in 
            easting_dem and northing_dem.
        contour_int     : int                       : Desired contour interval spacing in meters. Default is 20 m.
    RETURNS
        Contour intervals plotted on specified axis..
    '''
    # create contour intervals
    min_elev = np.nanmin(elev_dem)
    max_elev = np.nanmax(elev_dem)
    contour_levels = np.arange(contour_int * np.floor(min_elev/contour_int),
                               contour_int * np.ceil(max_elev/contour_int) + contour_int, 
                               step=contour_int)

    # add contours to figure
    ax.contour(easting_dem, northing_dem, elev_dem, 
               levels=contour_levels, 
               colors='grey', linewidths=0.8)

    return
