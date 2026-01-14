#!/usr/bin/python3
'''Utility functions for loading data.'''

import glob, obspy, os, glob
import pandas as pd
import numpy as np
import datetime
from obspy.core.util import AttribDict
from obspy.core.utcdatetime import UTCDateTime
#from obspy.geodetics.base import gps2dist_azimuth

import settings
settings.set_paths('laptop')



def load_mseed_data(path_data, path_coords, array_str=None, 
                    time_start=None, time_stop=None, freq_min=None, freq_max=None, 
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
        time_start : datetime object
        time_stop : "YYYY-MM-DD" 

        filter_type : str : Optional. Obspy filter type. Includes 'bandpass', 
           highpass', and 'lowpass'.
        filter_options : dict : Optional. Obspy filter arguments. For 'bandpass', 
            contains freq_min and freq_max. For low/high pass, contains freq.
    RETURNS
        data : obspy stream : Stream of data traces for full array, or specified 
            Gems. Stats include assigned coordinates.
    '''
    # (1) read in the data
    # get list of all dates to load mseed files (saved by date)
    dates = pd.date_range(start=time_start.date(), end=time_stop.date(), freq='D')
    data = obspy.Stream()
    for date in dates:    
        date_str = date.strftime('%Y-%m-%d')
        if array_str == None:   # if no array is specified, get all arrays
            path_mseed = os.path.join(path_data, f"{date_str}*.mseed")
        else:   # just get specified array
            path_mseed = os.path.join(path_data, f"{date_str}*.{array_str}.*.mseed")

        data += obspy.read(path_mseed)


    # (2) filter and detrend data if desired (needs to happen before data.merge)
    if freq_min != None and freq_max != None: # bandpass
        #data = data.detrend('linear')
        data = data.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
    elif freq_min != None:   # high pass
        data = data.detrend('linear')
        data = data.filter('highpass', freq=freq_min, corners=6)
    elif freq_max != None:   # low pass
        #data = data.detrend('linear')
        data = data.filter('lowpass', freq=freq_max)

    # (3) merge dates and trim data
    data = data.merge(method=0) # merge, discarding overlaps and leaving gaps
    data = data.trim(starttime=UTCDateTime(time_start), endtime=UTCDateTime(time_stop), keep_empty_traces=False)
    
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
        gem_include_station = coords[coords['SN'].isin(gem_include)]['Station'].to_list()
        data = obspy.Stream([trace for trace in data.traces if trace.stats['station'] in gem_include_station])
    if gem_exclude != None:
        # convert gem SNs (in exclude list) to stations
        gem_exclude_station = coords[coords['SN'].isin(gem_exclude)]['Station'].to_list()
        data = obspy.Stream([trace for trace in data.traces if trace.stats['station'] not in gem_exclude_station])

    return data



def filename_beamform(array_str, date, freq_min, freq_max):
    freq_str = f"{freq_min}-{freq_max}Hz"
    date_str = date.strftime('%Y-%m-%d')
    file_str = f"{array_str}_{date_str}_{freq_str}"
    path_processed = os.path.join(settings.path_processed, f"processed_output_{file_str}.pkl")
    return file_str, path_processed


def load_beamform_output(array_str, time_start, time_stop, freq_min, freq_max, slow_min=1/0.40, slow_max=1/0.30):
    '''
    Load in beamforming results for an array and frequency band in a specified time range. Remove points with 
    unreasonable slownesses, if desired.
    INPUTS
    RETURNS
    '''
    # load mulitple files if needed (otherwise, will just load 1)
    beamform_output = []
    dates = pd.date_range(start=time_start.date(), end=time_stop.date(), freq='D')
    for date in dates:    
        _, path_processed = filename_beamform(array_str, date, freq_min, freq_max)
        outputi = pd.read_pickle(path_processed)
        beamform_output.append(outputi)
    beamform_output = pd.concat(beamform_output)

    # extract times of interest
    beamform_output = beamform_output.loc[time_start : time_stop]

    # constrain data to only points with reasonable slownesses
    if slow_min != None and slow_max != None:
        beamform_output = beamform_output[beamform_output["Slowness"].between(slow_min, slow_max)]

    return beamform_output