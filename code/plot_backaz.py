
import datetime, pytz, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

import settings
settings.set_paths('laptop')
import utils.load_data


def plot_single_array_backaz():
    array_str = 'SC'
    time_start = datetime.datetime(2025, 10, 6, 0, 0, 0, tzinfo=pytz.timezone('UTC'))
    time_stop = datetime.datetime(2025, 10, 11, 23, 59, 59, tzinfo=pytz.timezone('UTC'))

    freq_min = 4
    freq_max = 8

    file_str, path_processed = utils.load_data.filename_beamform(array_str, time_start, freq_min, freq_max)

    beamform_output = utils.load_data.load_beamform_output(array_str, time_start, time_stop, 
                                                    freq_min=freq_min, freq_max=freq_max, 
                                                    slow_min=2.5, slow_max=3.5)


    plot_backaz(beamform_output, settings.path_home, 
                            f"{array_str} Array, Filtered {freq_min} to {freq_max} Hz", 
                            file_str)
    
    return



def plot_beamform_results():

    array_list = ['NW', 'NC', 'NE', 'SC']
    freq_list = [(4, 8)]
    
    time_start = datetime.datetime(2025, 10, 10, 16, 0, 0, tzinfo=pytz.timezone('UTC'))
    time_stop = datetime.datetime(2025, 10, 11, 5, 0, 0, tzinfo=pytz.timezone('UTC'))

    date = time_start.strftime("%Y-%m-%d")
    array_elements = utils.load_data.num_array_elements()

    for freq_min, freq_max in freq_list:
        print(f'{freq_min}-{freq_max} Hz')

        for array_str in array_list:
            print(array_str)
            beamform_output = utils.load_data.load_beamform_output(array_str, time_start, time_stop, 
                                                            freq_min=freq_min, freq_max=freq_max, 
                                                            slow_min=0, slow_max=5)
            beamform_output_filtered = utils.load_data.load_beamform_output(array_str, time_start, time_stop, 
                                                            freq_min, freq_max, 
                                                            slow_min=2.5, slow_max=3.5)
            
            font_label = 16
            font_title = 20
            alpha = 0.2

            fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=[12, 8])

            # plot backaz
            ax[0].scatter(beamform_output.index, beamform_output['Backaz'],
                          color='black', alpha=alpha)
            ax[0].scatter(beamform_output_filtered.index, beamform_output_filtered['Backaz'],
                          color='black')
            ax[0].set_ylabel("Backazimuth ($^o$)", fontsize=font_label)
            ax[0].set_ylim([0, 360])
            ax[0].set_yticks(ticks=np.arange(0, 360+90, 90))
            ax[0].tick_params(labelsize=font_label)

            # plot slowness
            ax[1].scatter(beamform_output.index, beamform_output['Slowness'],
                          color='black', alpha=alpha)
            ax[1].scatter(beamform_output_filtered.index, beamform_output_filtered['Slowness'],
                          color='black')
            ax[1].set_ylabel("Slowness (km/s)\n", fontsize=font_label)
            ax[1].set_ylim([0.5, 5])
            ax[1].set_yticks(ticks=np.arange(1, 6, step=1))
            ax[1].tick_params(labelsize=font_label)

            # calculate adjusted semblance
            N = array_elements[array_str][date]
            adj_semblance = lambda output: (N / (N-1)) * (output['Semblance'] - 1/N)

            # plot adjusted semblance
            ax[2].scatter(beamform_output.index, adj_semblance(beamform_output),
                          color='black', alpha=alpha)
            ax[2].scatter(beamform_output_filtered.index, adj_semblance(beamform_output_filtered),
                          color='black')
            ax[2].set_ylabel("Adjusted\nSemblance", fontsize=font_label)
            ax[2].set_yscale('log')
            ax[2].set_ylim([0.005, 1])
            ax[2].set_yticks(ticks=[0.01, 0.1, 1], labels=[0.01, 0.1, 1])
            ax[2].tick_params(labelsize=font_label)

            # format x-axis
            ax[2].set_xlabel(f"Local Time on {date}", fontsize=font_label)
            tick_spacing = 1
            ax[2].xaxis.set_major_locator(mdates.HourLocator(byhour=range(24), interval=tick_spacing, tz=pytz.timezone("US/Pacific")))
            ax[2].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=pytz.timezone("US/Pacific")))
            fig.autofmt_xdate()

            # format plot
            fig.suptitle(f'{array_str}', fontsize=font_title)
            fig.tight_layout()
            fig.savefig(os.path.join(settings.path_figures, "beamform_results",
                                     f"{array_str}_{date}_{freq_min}-{freq_max}Hz.png"),
                                     dpi=500)
            plt.close()

    return





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
    if fig == None and ax == None:
        fig, ax = plt.subplots(1, 1, figsize=[24, 6], tight_layout=True)

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
    plt.tight_layout()

    if file_str == None:
        return fig, ax
    else: 
        # save figure
        plt.savefig(os.path.join(path_home, "figures", f"backaz_{file_str}.png"), dpi=500)
        plt.close()
        return


if __name__ == '__main__':
    #plot_single_array_backaz()

    plot_beamform_results()