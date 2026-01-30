#!/usr/bin/python3

import os, datetime
import numpy as np
import pandas as pd
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal.array_analysis import array_processing
#from concurrent.futures import ProcessPoolExecutor
#from obspy.core.util import AttribDict
from matplotlib.dates import num2date
import pytz

import settings
import utils.load_data as load_data
#import utils.plot as plot

settings.set_paths('laptop')

def main(path_home, path_coords, array_str,
         gem_include=None, gem_exclude=None,
         time_start=None, time_stop=None, 
         freq_min=None, freq_max=None):

    path_data = os.path.join(path_home, "data", "mseed")
    file_str, path_processed = load_data.filename_beamform(array_str, time_start, freq_min, freq_max)

    # print progress to log file
    print(f"{datetime.datetime.now()} \t\t Loading and Filtering Data ")
    print("    "+path_processed)

    # load data
    data = load_data.load_mseed_data(path_data, path_coords=path_coords,
                     array_str=array_str,
                           gem_include=gem_include, gem_exclude=gem_exclude, 
                           time_start=time_start, time_stop=time_stop,
                           freq_min=freq_min, freq_max=freq_max)
    
    # print progress to log file
    #with open(os.path.join(path_home, "code", "log", "pylog.txt"), "a") as f:
    print(f"{datetime.datetime.now()} \t\t Processing Data ")
    print("    "+path_processed)#, file=f)

    # fiter and beamform 
    output = process_data(data, path_processed, 
                            #time_start=time_start, time_stop=time_stop, 
                            freq_min=freq_min, freq_max=freq_max)
    
    # print progress to log file
    #with open(os.path.join(path_home, "code", "log", "pylog.txt"), "a") as f:

    ## plot backaz time series
    #plot.plot_backaz(output, path_home, 
    #                        f"{array_str} Array, Filtered {freq_min} to {freq_max} Hz", file_str)
    
    return


def process_data(data, path_processed=None, time_start=None, time_stop=None, 
                 freq_min=None, freq_max=None):
    '''
    Run obspy array_processing() function to beamform data. Save in .npy format to specified 
    location. Returns output as np array, with backazimuths from 0-360.
    INPUTS
        data : obspy stream : merged data files
        path_processed : str : path and filename to save output as .npy
        time_start : obspy UTCDateTime : if specified, time to start beamforming. If not specified, 
            will use max start time from all Gems.
        time_stop : obspy UTCDateTime : if specified, time to end beamforming. If not specified, 
            will use min end time from all Gems.
    RETURNS
        output : pd dataframe : array with 5 rows of output from array_processing. Rows are: timestamp (UTC), 
            relative power (semblance), absolute power, backazimuth (from 0-360), and slowness (in s/km).
    '''
    #FIXME doc string above
    # if times are not provided, use max/min start and end times from gems
        # NOTE this can lead to data gaps if one gem started late!!
    if time_start == None:
        time_start = max([trace.stats.starttime for trace in data])
    else:
        time_start = UTCDateTime(time_start)
    if time_stop == None:
        time_stop = min([trace.stats.endtime for trace in data])
    else:
        time_stop = UTCDateTime(time_stop)

    
    #FIXME can clean this up when these change
    process_kwargs = dict(
        # slowness grid (in [s/km])
        sll_x=-4.0, slm_x=4.0, sll_y=-4.0, slm_y=4.0, sl_s=0.1,
        # sliding window
        win_len=60, win_frac=0.50,
        # frequency
        frqlow=freq_min, frqhigh=freq_max, prewhiten=0,
        # output restrictions
        semb_thres=-1e9, vel_thres=-1e9, timestamp='mlabday',
        stime=time_start, etime=time_stop)
    
    output = array_processing(stream=data, **process_kwargs)

    # correct backaz from 0 to 360 (instead of -180 to +180)
    output[:,3] = [output[i][3] if output[i][3]>=0 else output[i][3]+360 
                    for i in range(output.shape[0])]

    # save output as dataframe
    output = pd.DataFrame(data=output, 
                          columns=["Time", "Semblance", "Abs Power", "Backaz", "Slowness"])
    
    # save time steps as datetime types
    output["Time"] = num2date(output["Time"])
    # set index to time
    output = output.set_index("Time")

    # save output to pickle
    if path_processed != None:
        output.to_pickle(path_processed)
    return output




if __name__ == "__main__":

    #array_str = 'NW'
    #gem_exclude = None
    #array_str = 'NC'
    #gem_exclude = ['NC20']
    #array_str = 'NE'
    #gem_exclude = None
    #array_str = 'SC'
    #gem_exclude = #CHECK NOTEBOOK 16 JAN# ['SC05', 'SC14', 'SC20']

    freq_min = 4
    freq_max = 8
    
    array_str = 'SC'
    gem_exclude = ['SC05', 'SC14', 'SC20', 'SC09', 'SC03']
    time_start = datetime.datetime(2025, 10, 7, 0, 0, 0, tzinfo=pytz.timezone('UTC'))
    time_stop = datetime.datetime(2025, 10, 8, 0, 0, 0, tzinfo=pytz.timezone('UTC'))

    main(settings.path_home, settings.path_coords, array_str,
        gem_include=None, gem_exclude=gem_exclude,
        time_start=time_start, time_stop=time_stop, 
        freq_min=freq_min, freq_max=freq_max)

    
    #ndays = 6
    #time_start = datetime.datetime(2025, 10, 6, 0, 0, 0, tzinfo=pytz.timezone('UTC'))
    #time_start_list = [time_start + datetime.timedelta(days=x) for x in range(ndays)]
    #time_stop_list = [t + datetime.timedelta(days=1) for t in time_start_list]

    #for i in range(nhours):
    #    time_start = time_start_list[i]
    #    time_stop = time_stop_list[i]
    
    #    main(settings.path_home, settings.path_coords, array_str,
    #        gem_include=None, gem_exclude=gem_exclude,
    #        time_start=time_start, time_stop=time_stop, 
    #        freq_min=freq_min, freq_max=freq_max)




