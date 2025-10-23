#!/usr/bin/python3

import glob, obspy, os, glob, xmltodict, utm
import pandas as pd
import numpy as np
import datetime
from obspy.core.util import AttribDict
from obspy.geodetics.base import gps2dist_azimuth
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import obspy.signal.array_analysis

    
def load_mseed_data(path_data, path_coords, array_str=None, 
                    time_start=None, time_stop=None, freqmin=None, freqmax=None, 
                    gem_include=None, gem_exclude=None):
    '''
    Loads in mseed files between times of interest as an obspy stream. Assigns 
    coordinates from csv. Filters the data if desired, and ignores any specified Gem SNs.

        Loads all miniseed files in a specified directory into an obspy stream. 
        Assigns coordinates to all traces. If specified, filters data. If specified, 
        only returns a subset of gems (otherwise, returns full array).
    INPUTS
        path_data   : str : Path to mseed data folder.
        path_coords : str : Path and filename of coordinates.
        gem_include : list of str : Optional. If specified, should list Gem station
            names to include in processing. Mutually exclusive with gem_exclude.
        gem_exclude : list of str : Optional. If specified, should list Gem station
            names to include in processing. Mutually exclusive with gem_include.
        time_start : 
        time_stop : "YYYY-MM-DD" 

        filter_type : str : Optional. Obspy filter type. Includes 'bandpass', 
            'highpass', and 'lowpass'.
        filter_options : dict : Optional. Obspy filter arguments. For 'bandpass', 
            contains freqmin and freqmax. For low/high pass, contains freq.
    RETURNS
        data : obspy stream : Stream of data traces for full array, or specified 
            Gems. Stats include assigned coordinates.
    '''
    # (1) read in the data
    #TODO make this have multiple dates and use glob.glob
    date_str = str(time_start.date)
    if array_str != None:
        path_mseed = os.path.join(path_data, f"*{date_str}*{array_str}*.mseed")
    else:
        path_mseed = os.path.join(path_data, f"*{date_str}*.mseed")

    data = obspy.read(path_mseed)

    # (2) filter data if desired (needs to happen before data.merge)
    if freqmin != None and freqmax != None: # bandpass
        data = data.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
    elif freqmin != None:   # high pass
        data = data.filter('highpass', freq=freqmin)
    elif freqmax != None:   # low pass
        data = data.filter('lowpass', freq=freqmax)

    # (3) merge dates and trim data
    data = data.merge(method=0) # merge, discarding overlaps and leaving gaps
    data = data.trim(starttime=time_start, endtime=time_stop, keep_empty_traces=False)
    
    # (4) add coordinates to traces
    coords = pd.read_csv(path_coords)
    coords["SN"] = coords["SN"].astype(str)     # SN of gem (make sure csv has correct header)
    # get rid of any stations that don't have coordinates
    data = obspy.Stream([trace for trace in data.traces 
                    if trace.stats['station'] in coords['Station'].to_list()])
    # assign coordinates to stations
    for _, row in coords.iterrows():
        for trace in data.select(station=row["Station"]):
            trace.stats.coordinates = AttribDict({
                'latitude': row["Latitude"],
                'longitude': row["Longitude"],
                'elevation': row["Elevation"] }) 
    
    # (5) filter by gem station ID (only use specified subset)
    if gem_include != None:
        data = obspy.Stream([trace for trace in data.traces if trace.stats['station'] in gem_include])
    if gem_exclude != None:
        data = obspy.Stream([trace for trace in data.traces if trace.stats['station'] not in gem_exclude])

    return data


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

    # constrain data to only plot points with slownesses near 3 s/km
    slow_min = 2.5
    slow_max = 3.5
    output = output[output["Slowness"].between(slow_min, slow_max)]

    # create figure
    if fig == None and ax ==None:
        fig, ax = plt.subplots(1, 1, figsize=[7, 5], tight_layout=True)

    im = ax.scatter(output.index, output['Backaz'], c=output["Semblance"],
                    alpha=0.7, edgecolors='none', cmap='plasma',
                    vmin=min(output["Semblance"]), vmax=max(output["Semblance"]))
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Semblance")

    # format y-axis
    ax.set_ylabel("Backazimuth [$^o$]")
    ax.set_ylim([0, 360])
    ax.set_yticks(ticks=np.arange(0, 360+60, 60))

    # format x-axis
    ax.set_xlabel("Local Time")
    ax.set_xlim([output.index.min(), output.index.max()])
    hours_num = (output.index.max() - output.index.min()).total_seconds() / 3600
    tick_spacing = 1#int(np.ceil((hours_num / 15))) # make x-axis look nice (good number of ticks)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(24), interval=tick_spacing))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M", tz="US/Pacific"))
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
