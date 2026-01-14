#!/usr/bin/python3
'''Utility functions for plotting.'''

import obspy, os
import numpy as np
#import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import obspy.signal.array_analysis


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



def plot_backaz(output, path_home, subtitle_str, file_str=None, fig=None, ax=None):
    #TODO change path_home to path_figures
    '''
    Plot backazimuth over time from output of array processing. 
    INPUTS:
        output : pandas df : Result from beamforming with the columns 
            Time (datetime), Semblance, Abs Power, Backaz (0-360), and Slowness.
        path_home : str : Path to main dir. Figure will be saved in "figures" subdir.
        subtitle_str : str : Subtitle for plot. Usually contains bandpass frequencies 
            (e.g. "Filtered 24-32 Hz")
        file_str : str or None : String to append on end of filename to uniquely save figure 
            (e.g. "24.0_32.0"). If None, function returns a handle to the figure and axes, and does 
            NOT save the figure. 
    RETURNS:
        If file_str=None, returns handle to the figure and axes. Figure is NOT saved.
        Otherwise, figure is saved as path_home/figures/backaz_{file_str}.png
    '''
    # sort by ascending semblance so brightest points are plotted on top
    output = output.sort_values(by="Semblance", ascending=True)


    # create figure
    if fig == None and ax ==None:
        fig, ax = plt.subplots(1, 1, figsize=[7, 5], tight_layout=True)

    im = ax.scatter(output.index, output['Backaz'], c=output["Semblance"],
                    alpha=0.7, edgecolors='none', cmap='plasma',
                    vmin=0, vmax=1)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(['0', '0.5', '1'])
    cbar.ax.minorticks_off()
    cbar.set_label("Semblance")

    # format y-axis
    ax.set_ylabel("Backazimuth [$^o$]")
    ax.set_ylim([0, 360])
    ax.set_yticks(ticks=np.arange(0, 360+60, 60))

    # format x-axis
    ax.set_xlabel("Time (UTC)")
    ax.set_xlim([output.index.min(), output.index.max()])
    hours_num = (output.index.max() - output.index.min()).total_seconds() / 3600
    tick_spacing = 2#int(np.ceil((hours_num / 15))) # make x-axis look nice (good number of ticks)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(24), interval=tick_spacing))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M", tz="UTC"))#tz="US/Pacific"))
    fig.autofmt_xdate()

    # add titles
    fig.suptitle(f"Backazimuth")
    ax.set_title(subtitle_str, fontsize=10)

    if file_str == None:
        return fig, ax
    else: 
        # save figure
        plt.savefig(os.path.join(path_home, "figures", f"backaz_{file_str}.png"), dpi=500)
        plt.close()
        return


