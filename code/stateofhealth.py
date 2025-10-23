#!/usr/bin/python3

import os, math, utm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import settings
settings.set_paths('laptop')

array_str = "B"

# get all gem SNs from coords
coords = pd.read_csv(settings.path_coords)

coords = coords[coords["Network"] == array_str]
coords = coords.reset_index()

n_stations = len(coords)
fig1, ax1 = plt.subplots(ncols=2, nrows=math.ceil(n_stations/2), sharex=True, sharey=True, figsize=[7, 14])
fig2, ax2 = plt.subplots(ncols=2, nrows=math.ceil(n_stations/2), sharex=True, sharey=True, figsize=[7, 14])

# set median x limits to ignore data not erased before deployment
xmin = []
xmax = []


for i, row in coords.iterrows():
    sn = row["SN"]
    axi1 = ax1.ravel()[i]
    axi2 = ax2.ravel()[i]
    axi1.set_title(sn)
    axi2.set_title(sn)

    try:
        soh = pd.read_csv(os.path.join(settings.path_home, "data", array_str, "metadata", f"{sn}metadata_000.txt"))
    except:
        # for stations that did not gemconvert properly
        continue

    axi1.plot(soh['t'], soh['batt'], 'b-')
    axi2.plot(soh['t'], soh['temp'],  'r-')

    xmin.append(soh['t'].min())
    xmax.append(soh['t'].max())

# set x lim
axi1.set_xlim([np.median(xmin), np.median(xmax)])
axi2.set_xlim([np.median(xmin), np.median(xmax)])

fig1.suptitle("Battery Voltage (V)")
fig2.suptitle("Temperature (C)")

fig1.tight_layout()
fig2.tight_layout()

fig1.savefig(os.path.join(settings.path_figures, "stateofhealth", f"{array_str}_battery_voltage.png"))
fig2.savefig(os.path.join(settings.path_figures, "stateofhealth", f"{array_str}_temperature.png"))

plt.close()
#plt.show()




