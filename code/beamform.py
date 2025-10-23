#!/usr/bin/python3

import os, datetime
import numpy as np
import pandas as pd
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal.array_analysis import array_processing
#from concurrent.futures import ProcessPoolExecutor
#from obspy.core.util import AttribDict
from matplotlib.dates import num2date

import utils, settings

def main(path_home, path_coords, array_str,
         gem_include=None, gem_exclude=None,
         time_start=None, time_stop=None, 
         freqmin=None, freqmax=None):

    path_data = os.path.join(path_home, "data", array_str, "mseed")
    filt_freq_str = f"{freqmin}_{freqmax}"
    filt_date_str = f"{time_start.year}-{time_start.month}-{time_start.day}"
    file_str = f"{array_str}_{filt_date_str}_{filt_freq_str}"

    path_processed = os.path.join(path_home, "data", array_str, "processed", 
                    f"processed_output_{file_str}.pkl")

    # print progress to log file
    #with open(os.path.join(path_home, "code", "log", "pylog.txt"), "a") as f:
    print((f"{datetime.datetime.now()} \t\t Loading and Filtering Data ")+
            (f"({file_str})"))
    print("    "+path_data)

    # load data
    data = utils.load_mseed_data(path_data, path_coords=path_coords,
                     array_str=array_str,
                           gem_include=gem_include, gem_exclude=gem_exclude, 
                           time_start=time_start, time_stop=time_stop,
                           freqmin=freqmin, freqmax=freqmax)
    
    # print progress to log file
    #with open(os.path.join(path_home, "code", "log", "pylog.txt"), "a") as f:
    print((f"{datetime.datetime.now()} \t\t Processing Data ")+
            (f"({file_str})"))#, file=f)
    print("    "+path_processed)#, file=f)

    # fiter and beamform 
    output = process_data(data, path_processed, 
                            time_start=None, time_stop=None, 
                            freqmin=freqmin, freqmax=freqmax)
    


    # print progress to log file
    #with open(os.path.join(path_home, "code", "log", "pylog.txt"), "a") as f:
    print((f"{datetime.datetime.now()} \t\t Plotting Backazimuth ")+ 
            (f"({file_str})"))
    print("    "+os.path.join(path_home, "figures", f"backaz_{file_str}.png"))

    # plot backaz time series
    utils.plot_backaz(output, path_home, 
                            f"{array_str} Array, Filtered {freqmin} to {freqmax} Hz", file_str)
    
    return


def process_data(data, path_processed=None, time_start=None, time_stop=None, freqmin=None, freqmax=None):
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
    if time_start == None:
        time_start = max([trace.stats.starttime for trace in data])
    if time_stop == None:
        time_stop = min([trace.stats.endtime for trace in data])

    
    #FIXME can clean this up when these change
    process_kwargs = dict(
        # slowness grid (in [s/km])
        sll_x=-4.0, slm_x=4.0, sll_y=-4.0, slm_y=4.0, sl_s=0.1,
        # sliding window
        win_len=60, win_frac=0.50,
        # frequency
        frqlow=freqmin, frqhigh=freqmax, prewhiten=0,
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

    settings.set_paths('laptop')
    array_str = 'A'
    gem_exclude = None


    time_start = UTCDateTime(2025, 10, 6, 0, 0, 0)
    time_stop = UTCDateTime(2025, 10, 11, 0, 0, 0)
    freqmin = 2
    freqmax = 10

    
    main(settings.path_home, settings.path_coords, array_str,
         gem_include=None, gem_exclude=gem_exclude,
         time_start=time_start, time_stop=time_stop, 
         freqmin=freqmin, freqmax=freqmax)




