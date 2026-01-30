#!/usr/bin/python3
'''Plot array geometry, array response, and network response.'''
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import pandas as pd
import numpy as np
import scipy as sci
import os, geopandas, itertools, math, datetime

import settings, utils.plot, utils.load_data
settings.set_paths('laptop')


#TODO MAKE SURE BAD STATIONS ARE REMOVED****


def array_response():
    # import coordinates
    coords = pd.read_csv(settings.path_coords)

    # initialize figure
    fig, ax = plt.subplots(ncols=4, nrows=2, figsize=[12,5])
    # figsize=[6,10] for vertical plot

    array_list = ['NW', 'NC', 'NE', 'SC']
    array_name = ['NW', 'NC', 'NE', 'SC']

    font_ax = 12
    font_title = 14

    # remove SN 245 (did not record)
    #coords.loc[coords['SN'] == 245, 'Network'] = 0

    for i, array_str in enumerate(array_list):
        # get coordinates for this array
        coords_array = coords[coords['Network'] == array_str]

        # plot array geometry
        ax[0,i].plot(coords_array['Easting'], coords_array['Northing'], 
                     'k^', markersize=10)

        # set tick labels relative to array center
        center = [coords_array['Easting'].mean(), coords_array['Northing'].mean()]
        tick_max = 100  # distance in m
        tick_num = 5    # number of x/y ticks
        ax[0,i].set_xticks(ticks=np.linspace(center[0]-tick_max, center[0]+tick_max, tick_num),
                        labels = np.linspace(-tick_max, tick_max, tick_num, dtype=int))
        ax[0,i].set_yticks(ticks=np.linspace(center[1]-tick_max, center[1]+tick_max, tick_num),
                        labels = np.linspace(-tick_max, tick_max, tick_num, dtype=int))

        # format axes
        ax[0,i].set_title(array_name[i], fontsize=font_title)
        #ax[0,0].set_title(r"$\bf{Array\ Geometry}$"+f"\n{array_name[0]}")
        ax[0,i].tick_params(labelsize=font_ax)
        ax[0,i].set_xlabel('Distance (m)', fontsize=font_ax)
        ax[0,i].set_ylabel('Distance (m)', fontsize=font_ax)
        ax[0,i].set_aspect('equal')
        ax[0,i].grid(True)

        # calculate and plot array response
        res = utils.plot.plot_array_response(coords_array, flim=15, ax=ax[1,i])

        # plot circles at relevant freqs
        plot_wavenum_circles(freqs=[4,8], ax=ax[1,i])

        # format axes
        ax[1,i].set_title(array_name[i], fontsize=font_title)
        #ax[1,0].set_title(r"$\bf{Array\ Response}$"+f"\n{array_name[1]}")
        ax[1,i].tick_params(labelsize=font_ax)
        ax[1,i].set_xlabel("$k_x$ (rad/km)", fontsize=font_ax)
        ax[1,i].set_ylabel("$k_y$ (rad/km)", fontsize=font_ax)
        ax[1,i].set_aspect('equal')
        ax[1,i].sharex(ax[1, i-1])
        ax[1,i].sharey(ax[1,i-1])
        ax[1,i].set_xticks([-250,0,250])
        ax[1,i].set_yticks([-250,0,250])

        # set colorbars
        cbar = fig.colorbar(res, ax=ax[1,i], 
                            fraction=0.046, pad=0.03)
        cbar.set_ticks([0.0, 1.0])
        cbar.ax.tick_params(labelsize=font_ax)
        cbar.set_label("Semblance", fontsize=font_ax)
    
    fig.tight_layout()
    # save figure
    fig.savefig(os.path.join(settings.path_figures, "fig_array_response_horiz.png"), dpi=500)
    return


def plot_wavenum_circles(freqs, ax):
    '''
    Plot wavenumber circles on array reponse plot.
    '''
    if type(freqs) is not list:
        freqs = [freqs]
    for f in freqs:
        r = 2*np.pi*f/0.343
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(r*np.sin(theta), r*np.cos(theta), '--', color='grey')
    return





def network_response(remove_array=None):
    
    # define bounds in UTM coords (min_east, max_east, min_north, max_north)
    bounds_region = [668500, 672200, 5199000, 5201300]
    #bounds_region = [668529, 673367, 5197963, 5201585]

    # set up figure
    fig, ax = plt.subplots(1, 1, figsize=[10,7])
    ax.set_xlim(bounds_region[0], bounds_region[1])
    ax.set_ylim(bounds_region[2], bounds_region[3])
    ax.ticklabel_format(style='plain')
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    # create scale bar
    scale_bar = 500 # meters
    ax.hlines(y=bounds_region[2]+100, xmin=bounds_region[0]+100, xmax=bounds_region[0]+100+scale_bar, 
              color='black', linewidth=4)
    ax.text(y=bounds_region[2]+130, x=bounds_region[0]+160, s=f'{scale_bar} m', fontsize=16,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])

    # read in DEM and plot contour lines
    easting_dem, northing_dem, elev_dem, crs = utils.plot.read_dem(settings.path_gis_dem, bounds_region)
    utils.plot.plot_contours_from_dem(ax, easting_dem, northing_dem, elev_dem, 
                                      contour_int=20)

    # read in shapefile and plot intended burn area
    gdf = geopandas.read_file(settings.path_gis_unit)
    gdf = gdf.to_crs(crs)
    gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=1.5)

    # create meshgrid of all locations for network response
    easting = np.arange(bounds_region[0], bounds_region[1], step=25)
    northing = np.arange(bounds_region[2], bounds_region[3], step=25)
    [xgrid, ygrid] = np.meshgrid(easting, northing)
    # get DEM elevations on same meshgrid
    interpolator = sci.interpolate.RegularGridInterpolator(points=(northing_dem[:,1], easting_dem[0,:]),
                                           values=elev_dem,
                                           method='linear',
                                           bounds_error=False, fill_value=np.nan)
    zgrid = interpolator(np.column_stack([ygrid.ravel(), xgrid.ravel()]))
    zgrid = zgrid.reshape(xgrid.shape)


    # read in station coordinates
    coords = pd.read_csv(settings.path_coords)
    #TODO REMOVE BAD SENSORS/ CHANGE N??

    # remove cameras (only plot infrasound)
    coords = coords[coords['Network'] != 'T']
    #NOTE remove one station at a time for sensitivity analysis
    if remove_array != None:
        coords = coords[coords['Network'] != remove_array]
        coords = coords.reset_index()

    for array_str in coords['Network'].unique():
        coords_sub = coords[coords['Network'] == array_str]
        ax.plot(coords_sub['Easting'], coords_sub['Northing'], 'k^', markersize=12)
        # label arrays
        ax.text(coords_sub['Easting'].max(), coords_sub['Northing'].max(), 
                array_str, fontsize=16, zorder=100, 
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])

    # calculate distances btw all stations and all grid points
    all_dists = calculate_grid_to_station_distances(coords, xgrid, ygrid, zgrid)

    # subsample meshgrid for plotting
    spacing = 10
    xgrid_sub = xgrid[0::spacing, 0::spacing]
    ygrid_sub = ygrid[0::spacing, 0::spacing]
    #ax.plot(xgrid_sub, ygrid_sub, '.', color='indigo')


    # loop through each point in the subgrid
    for k, (xtmp, ytmp) in enumerate(zip(xgrid_sub.ravel(), ygrid_sub.ravel())):
        print(f"k={k} of {np.size(xgrid_sub)}")

        # find index of grid point in the full meshgrid
        idxx = np.where(np.isin(easting, xtmp))
        idxy = np.where(np.isin(northing, ytmp))
        
        # get distances from current reference source to all sensors
        ref_dist = np.squeeze(all_dists[idxy, idxx, :])

        # calc residuals between reference source and all other possible sources
        #TODO make this a function (with networks removable)
        total_res = 0
        for array_str in coords['Network'].unique():
            res = calculate_station_residuals(array_str, coords, ref_dist, all_dists)
        
            # divide by number of elements
            #date = datetime.datetime(2025, 10, 8)
            #N = utils.load_data.num_array_elements()[array_str][date]
            if array_str == 'NW' or 'NE':
                N = 14
            else:
                N = 25
            res = res / N
            total_res += res

        # calc RMS
        #total_res = (res_nw/14 + 
        #             res_nc/25 +  
        #             res_ne/14 +  
        #             res_sc/25)  
        dist_res = np.sqrt(total_res)
        # get time of arrival residuals (in sec)
        time_res = dist_res/343
        #print(np.nanmin(time_res), np.nanmax(time_res))

        if remove_array == None:
            c1 = 'darkgreen'
        else:
            c1 = 'darkblue'
        ax.contour(easting, northing, time_res, 
                   levels=np.arange(0.01, 0.03, 0.01), 
                   colors=[c1], alpha=0.5)
        ax.contour(easting, northing, time_res, 
                   levels=np.linspace(0.01, 0.11, 1), 
                   colors=[c1])

    ax.plot(0,0,'-', color=c1, linewidth=3,
               label='0.01 s Time Residual')
    ax.plot(1,1,'-', color=c1, linewidth=3, alpha=0.5,
                      label='0.02 s Time Residual')
    handles, _ = plt.gca().get_legend_handles_labels()
    ax.legend(handles=handles, fontsize=14, loc='lower right',
                  framealpha=1)
    
    fig.tight_layout()

    if remove_array == None:
        fig.suptitle("Wenas Network Response", fontsize=16)
        fig.savefig(os.path.join(settings.path_figures, "fig_network_response.png"), dpi=500)
    else:
        fig.suptitle(f"Wenas Network Response, {remove_array} Removed", fontsize=16)
        fig.savefig(os.path.join(settings.path_figures, 
                                 f"fig_network_response_{remove_array}_removed.png"), dpi=500)
    #plt.show()
    return


def calculate_grid_to_station_distances(coords, xgrid, ygrid, zgrid):
    '''
    Calculate the Euclidean distance from every point on a meshgrid to each station.
    INPUTS
        coords  : pd dataframe      : Contains columns Easting, Northing, and Elevation
        xgrid   : np array (nx, ny) : Meshgrid of all easting coordinates in region of interest.
        ygrid   : np array (nx, ny) : Meshgrid of all northing coordinates in region of interest.
        zgrid   : np array (nx, ny) : Elevations at all easting and northing coordinates.
    RETURNS
        dists   : np array (nx, ny, nstations)  : Array of Euclidean distance between 
            each sensor and each point on the meshgrid. Each slice [:,:,k] gives the 
            distance from sensor k to eahc point on the meshgrid.
    '''
    # get grid points in correct shape to manipulate
    xg = xgrid[:,:,np.newaxis]
    yg = ygrid[:,:,np.newaxis]
    zg = zgrid[:,:,np.newaxis]

    # get station coordinates
    xstation = coords["Easting"].to_numpy()[np.newaxis,np.newaxis,:]
    ystation = coords["Northing"].to_numpy()[np.newaxis,np.newaxis,:]
    zstation = coords["Elevation"].to_numpy()[np.newaxis,np.newaxis,:]

    # calculate distance between grid points and station coordinates
    all_dists = np.sqrt((xg - xstation)**2 + (yg - ystation)**2 + (zg - zstation)**2)
    return all_dists



def calculate_station_residuals(array_str, coords, ref_dist, all_dists):
    '''
    Calculate distance residuals between a reference source and all other possible sources.
    INPUTS
        ref_dist    : 
        all_dists   : 
    RETURNS
    '''
    idx = coords.index[coords['Network']==array_str]
    idx_start = idx[0]
    idx_stop = idx[-1]
    # make list of indices of all stations in this array
    idxs = list(range(idx_start, idx_stop+1))
    # get all pairs of stations in this array
    station_pairs = list(itertools.combinations(idxs, 2))
    # create empty array for distance residuals
    residuals = np.zeros(shape=(all_dists.shape[0], all_dists.shape[1]))

    # loop through all station pairs in the array
    for i, j in station_pairs:
        # get difference in dist from current reference source
        ref_diff = ref_dist[i] - ref_dist[j]
        # get difference in dist from all other sources
        all_diff = all_dists[:,:,i] - all_dists[:,:,j]
        # calculate residual and store in array
        residuals += (ref_diff - all_diff)**2
    return residuals



if __name__ == "__main__":
    array_response()

    #network_response(remove_array=None)