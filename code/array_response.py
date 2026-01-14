#!/usr/bin/python3
'''Plot array geometry and array response.'''
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

import utils.plot as plot, settings
settings.set_paths('laptop')


def main():
    # import coordinates
    coords = pd.read_csv(settings.path_coords)

    # initialize figure
    fig, ax = plt.subplots(ncols=2, nrows=4, figsize=[7,12])

    array_list = ['A', 'B', 'C', 'D']
    array_name = ['Array A', 'Array B', 'Array C', 'Array D']

    # remove SN 245 (did not record)
    #coords.loc[coords['SN'] == 245, 'Network'] = 0

    for i, array_str in enumerate(array_list):
        coords_array = coords[coords['Network'] == array_str]

        # plot array geometry
        ax[i,0].plot(coords_array['Easting'], coords_array['Northing'], 'k^', markersize=10)

        # set tick labels relative to array center
        center = [coords_array['Easting'].mean(), coords_array['Northing'].mean()]
        tick_max = 100
        tick_num = 5
        ax[i,0].set_xticks(ticks=np.linspace(center[0]-tick_max, center[0]+tick_max, tick_num),
                        labels = np.linspace(-tick_max, tick_max, tick_num, dtype=int))
        ax[i,0].set_yticks(ticks=np.linspace(center[1]-tick_max, center[1]+tick_max, tick_num),
                        labels = np.linspace(-tick_max, tick_max, tick_num, dtype=int))

        # label sensors with just station number
        for _, row in coords_array.iterrows():
            ax[i,0].annotate(text=f'{row["Station"].replace(array_str, "")}',
                            xy=(row['Easting'], row['Northing']),
                            xytext=(5, -5), textcoords='offset points')
        # format axes
        ax[i,0].set_title(array_name[i])
        ax[0,0].set_title(r"$\bf{Array\ Geometry}$"+f"\n{array_name[0]}")
        ax[i,0].set_xlabel('Distance (m)')
        ax[i,0].set_ylabel('Distance (m)')
        ax[i,0].set_aspect('equal')
        ax[i,0].grid(True)


        # calculate and plot array response
        res = plot.plot_array_response(coords_array, flim=20, ax=ax[i,1])

        # plot circles at relevant freqs
        plot_wavenum_circles(freqs=[2,10,20], ax=ax[i,1])

        # format axes
        ax[i,1].set_title(array_name[i])
        ax[i,1].set_title(r"$\bf{Array\ Response}$"+f"\n{array_name[1]}")
        ax[i,1].set_xlabel("$k_x$ (rad/km)")
        ax[i,1].set_ylabel("$k_y$ (rad/km)")
        ax[i,1].set_aspect('equal')
        ax[i,1].sharex(ax[i-1,1])
        ax[i,1].sharey(ax[i-1,1])
        ax[i,1].set_xticks([-250,0,250])
        ax[i,1].set_yticks([-250,0,250])

        # set colorbars
        cbar = fig.colorbar(res, ax=ax[i,1], fraction=0.046, pad=0.03)
        cbar.set_ticks([0.0, 0.5, 1.0])
        cbar.set_label("Semblance")

    fig.tight_layout()
    # save figure
    fig.savefig(os.path.join(settings.path_figures, "fig_array_response.png"), dpi=500)
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


if __name__ == "__main__":
    main()