#!/usr/bin/python3

import os, datetime
import numpy as np
import pandas as pd
import scipy as sci
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime, pytz
from obspy.signal.array_analysis import get_geometry, get_timeshift
from obspy.core.utcdatetime import UTCDateTime
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker


import settings
import utils.load_data as load_data
settings.set_paths('laptop')

print(datetime.datetime.now())

# parameters
time_start = datetime.datetime(2025, 10, 6, 20, 0, 0, tzinfo=pytz.timezone('UTC'))
#time_start = datetime.datetime(2025, 10, 7, 0, 0, 0, tzinfo=pytz.timezone('UTC'))
time_stop = datetime.datetime(2025, 10, 7, 4, 0, 0, tzinfo=pytz.timezone('UTC'))
array_str = 'A'
gem_exclude = None    
#array_str = 'B'
#gem_exclude = ['123']
#array_str = 'C'
#gem_exclude = None
#array_str = 'D'
#gem_exclude = ['005', '086', '368']



# import station coordinates and metadata
coords = pd.read_csv(settings.path_coords)

# load in waveforms
stream = load_data.load_mseed_data(settings.path_mseed, path_coords=settings.path_coords,
                    array_str=array_str,
                    gem_include=None, gem_exclude=gem_exclude, 
                    time_start=time_start, time_stop=time_stop,
                    freq_min=0.25, freq_max=None)   # detrend and highpass filter above 0.05 Hz (from Jake paper)
# TODO try 4 sec high pass



##spectra = np.empty([301, 9599, len(stream)])
##data = np.empty([2880001, len(stream)])
##for n, tr in enumerate(stream):
##    f, t, spectrum = sci.signal.spectrogram(x=tr.data, fs=100, 
##                                            nperseg=600, noverlap=300, 
##                                            detrend='constant', scaling='spectrum')
##    spectra[:,:,n] = spectrum
##    data[:,n] = tr.data
##
##
##spectra_mean = np.mean(spectra, axis=2)
##data_mean = np.mean(data, axis=1)
##
##fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 3]})
### plot mean
##ax[0].plot(stream[0].times(), data_mean, '-')
##
##
##im = ax[1].imshow(spectrum, 
##               cmap='inferno_r', 
##               norm=LogNorm(vmin=np.nanmin(spectra_mean), vmax=np.nanmax(spectra_mean)),
##               #vmin=np.nanmin(spectrum), 
##               #vmax=np.nanmax(spectrum),
##               origin='lower',
##               extent=[t.min(), t.max(), f.min(), f.max(),],
##               aspect='auto')
##plt.show()
##print(spectrum)
##plt.close()









################################################
# load in beamforming results
beamform_output = load_data.load_beamform_output(array_str, time_start, time_stop, 
                                                 freq_min=0.5, freq_max=2, slow_min=None, slow_max=None)


f_samp = stream[0].stats['sampling_rate']
geom = get_geometry(stream, coordsys='lonlat')#, return_center=True)
# return center and add center point to all elevations (to get no negative elevs)
#TODO maybe get array geom with east/northing coords from coord df instead
#( *** check this func and plot)





Y_stack_matrix = np.empty([3001, len(beamform_output)-1])
semblance_matrix = np.empty([3001, len(beamform_output)-1])
adj_semblance_matrix = np.empty([3001, len(beamform_output)-1])

for i, (t, row) in enumerate(beamform_output.iterrows()):
    print(i, datetime.datetime.now())
    if i+1 >= semblance_matrix.shape[1]:
        continue    #TODO hack to not deal with t, t[i+1] below for last entry
    
    # calculate slowness
    sx = row['Slowness'] * np.sin(np.deg2rad(row['Backaz']))
    sy = row['Slowness'] * np.cos(np.deg2rad(row['Backaz']))
    sz = row['Slowness'] * np.cos(np.deg2rad(90))   # TODO assuming indidence angle of 90 (horizontal)
    s = np.array([sx, sy, sz])
    #TODO make sure lat,lon vs lon,lat and might need to change order here

    # calculate time shifts for each station (eq 1)
    time_shifts = np.dot(geom, s)   # geom in km, s in s/km

    # get 30 s window of data corresponding with beamform output
    spectra = []
    for n, tr in enumerate(stream):
        # stream data already detrended and high-pass filter (over 0.05 Hz)
        slice = tr.slice(UTCDateTime(t), UTCDateTime(beamform_output.index[i+2]))

        try:
            slice = slice.detrend('linear')
        except:
            pass

        
        #TODO detrend slice
        
        
        #TODO change to 60 s window and step by 30 s still  
        npts = len(slice.data)
        data = slice.data * np.hamming(npts)      # window data

        # calculate fft
        freqs = sci.fft.fftfreq(npts, d=1/f_samp)
        freqs = freqs[:int(len(freqs)/2+1)]
        fft = sci.fft.rfft(data)
        # apply time shift (eq 2 innards)
        spectrum = fft * np.exp(2j * np.pi * freqs * time_shifts[n])
        # TODO neg or no negative 
        
        spectra.append(spectrum)
    
    # calculate stacked spectrum (eq 2 outards)
    spectra = np.array(spectra)
    Y_stack = np.mean(spectra, axis=0)

    # calculate semblance spectrum (eq 3)
    semblance = np.abs(Y_stack)**2 / np.mean(np.abs(spectra)**2, axis=0)

    # calculate adjusted semblance (eq 4)
    N = len(stream) # number of stations
    adj_semblance = (N / (N-1)) * (semblance - 1/N)

    if semblance.shape[0] != semblance_matrix.shape[0]:
        print(i)
        Y_stack_matrix[:,i] = np.zeros(3001)
        semblance_matrix[:,i] = np.zeros(3001)
        adj_semblance_matrix[:,i] = np.zeros(3001)
        continue    # TODO hack if there are time gaps in beamform_output, skip and leave this row blank

    Y_stack_matrix[:,i] = Y_stack
    semblance_matrix[:,i] = semblance
    adj_semblance_matrix[:,i] = adj_semblance

        

fig, ax = plt.subplots(1, 1)

times = beamform_output.index.to_numpy()
x = mdates.date2num(times)

# see what Y_stack**2 (D,E,F) looks like 
im = ax.imshow(adj_semblance_matrix, 
               cmap='inferno_r', 
               norm=LogNorm(vmin=0.01, vmax=1),
               #vmin=0, 
               #vmax=1,
               origin='lower',
               extent=[x.min(), x.max(), freqs.min(), freqs.max(),],
               aspect='auto')
cbar = fig.colorbar(im, label='Semblance', 
                    orientation='vertical', pad=0.02)
cbar.set_ticks([0.01, 0.3, 1])
cbar.set_ticklabels(['0.01', '0.3', '1'])
cbar.ax.minorticks_off()



ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
fig.autofmt_xdate()

ax.set_ylim([0, 20])


ax.set_ylabel('Frequency (Hz)')
fig.suptitle(f'Array {array_str} Semblance Spectrogram')

plt.tight_layout()
plt.savefig(os.path.join(settings.path_figures, 'test', 'test.png'), dpi=500)
#plt.show()











