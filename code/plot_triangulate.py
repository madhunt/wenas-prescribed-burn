#!/usr/bin/python3


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe
import geopandas, os, datetime, pytz
import matplotlib.dates as mdates

import settings, utils.plot
settings.set_paths('laptop')

def plot_triangualted_points(time_start, time_stop):


    # TODO change this so that its pulling from the full day file not a new file

    format = '%Y-%m-%dT%H:%M'
    filename = f'source_locs_{time_start.strftime(format)}-{time_stop.strftime(format)}_4-8Hz' 
    source_loc = pd.read_pickle(os.path.join(settings.path_processed, "triang_source_locations", f"{filename}.pkl"))


    #TODO get all of the dem/burn outline into a function in plot_utils
    bounds_region = [668500, 673000, 5199000, 5201300]


    # set up figure
    fig, ax = plt.subplots(1, 1, figsize=[10,6])
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

    # read in station coordinates
    #TODO make this an option to remove bad coords
    coords = pd.read_csv(settings.path_coords)
    # remove cameras (only plot infrasound)
    coords = coords[coords['Network'] != 'T']
    ax.plot(coords['Easting'], coords['Northing'], 'k^', markersize=12)
    # label arrays
    #TODO make this neat please
    coords_sub = lambda array_str: coords[coords['Network'] == array_str]
    ax.text(coords_sub('NW')['Easting'].max(), coords_sub('NW')['Northing'].max(),
            'NW', fontsize=16, zorder=100,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    ax.text(coords_sub('NC')['Easting'].max(), coords_sub('NC')['Northing'].max(),
            'NC', fontsize=16, zorder=100,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    ax.text(coords_sub('NE')['Easting'].max(), coords_sub('NE')['Northing'].max(),
            'NE', fontsize=16, zorder=100,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    ax.text(coords_sub('SC')['Easting'].max(), coords_sub('SC')['Northing'].max(),
            'SC', fontsize=16, zorder=100,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])

    #TODO above here is in setup image function


    #TODO check all this nonsense
    colors = plt.cm.Oranges(np.linspace(0, 1, len(source_loc)))
    norm = mpl.colors.Normalize(vmin=0, vmax=len(source_loc)-1)
    sm = mpl.cm.ScalarMappable(cmap=plt.cm.Oranges, norm=norm)
    sm.set_array([])
    # Add colorbar to the figure
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', 
                        shrink=0.8, aspect=30, pad=0.05)#, location='top')
    # Set colorbar tick positions and labels
    tick_positions = np.linspace(0, len(source_loc)-1, 5)
    
    # change tick labels to local time
    ticks = pd.date_range(time_start.astimezone(pytz.timezone('US/Pacific')), 
                          time_stop.astimezone(pytz.timezone('US/Pacific')), periods=5)
    tick_labels = ticks.strftime('%H:%M')
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels, fontsize=16)
    cbar.set_label(f'Local Time (UTC-7:00) on {time_start.strftime("%Y-%m-%d")}', rotation=0, fontsize=16)
    
    for i, (t, row) in enumerate(source_loc.iterrows()):
        print(f'{i}/{len(source_loc)}')
        # in a for loop for a potential future animation
        ax.plot(row['Easting'], row['Northing'],
                'o', markersize=8,
                markerfacecolor=colors[i], markeredgecolor='black')
        
    plt.tight_layout()
    #plt.show()
    plt.savefig(os.path.join(settings.path_figures, 'triang_source_locations', f'{filename}.png'))    

    return


def plot_array_contributions(time_start, time_stop):

    format = '%Y-%m-%dT%H:%M'
    filename = f'source_locs_{time_start.strftime(format)}-{time_stop.strftime(format)}_4-8Hz' 
    source_loc = pd.read_pickle(os.path.join(settings.path_processed, "triang_source_locations", f"{filename}.pkl"))

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, sharey=True, figsize=[10,8])

    date = time_start.strftime("%Y-%m-%d")
    font_label = 16
    font_title = 20

    array_list = ['NW', 'NC', 'NE', 'SC']
    for i, array_str in enumerate(array_list):
        ax[i].plot(source_loc.index, source_loc[f'{array_str} Ints'], 'ko')
        ax[i].set_title(f"{array_str}", fontsize=font_title)
        ax[i].set_ylabel("Number of\nIntersections", fontsize=font_label)
        ax[i].set_yticks(ticks=[0,1,2,3])
        ax[i].tick_params(labelsize=font_label)

        # calculate rolling average and plot
        win = 2*15 # rolling avg of 15 min
        weights = np.ones(win) / win
        moving_avg = np.convolve(source_loc[f'{array_str} Ints'].to_numpy(), weights, 'same')
        ax[i].plot(source_loc.index, moving_avg, 'r-')


    # format x-axis
    ax[3].set_xlabel(f"Local Time on {date}", fontsize=font_label)
    tick_spacing = 1
    ax[3].xaxis.set_major_locator(mdates.HourLocator(byhour=range(24), interval=tick_spacing, tz=pytz.timezone("US/Pacific")))
    ax[3].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=pytz.timezone("US/Pacific")))
    fig.autofmt_xdate()

    plt.tight_layout()
    fig_filename = f'int_contributions_{time_start.strftime(format)}-{time_stop.strftime(format)}_4-8Hz' 
    plt.savefig(os.path.join(settings.path_figures, 'triang_source_locations', f'{fig_filename}.png'))    



    return



if __name__ == "__main__":

    time_start = datetime.datetime(2025, 10, 10, 16, 0, 0, tzinfo=pytz.timezone('UTC'))
    time_stop = datetime.datetime(2025, 10, 11, 5, 0, 0, tzinfo=pytz.timezone('UTC'))

    plot_triangualted_points(time_start, time_stop)
    plot_array_contributions(time_start, time_stop)