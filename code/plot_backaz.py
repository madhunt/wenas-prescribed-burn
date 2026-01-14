
import datetime, pytz

import settings
settings.set_paths('laptop')
from utils.plot import plot_backaz
from utils.load_data import load_beamform_output, filename_beamform


array_str = 'A'
time_start = datetime.datetime(2025, 10, 6, 0, 0, 0, tzinfo=pytz.timezone('UTC'))
time_stop = datetime.datetime(2025, 10, 7, 0, 0, 0, tzinfo=pytz.timezone('UTC'))

freq_min = 2
freq_max = 4



file_str, path_processed = filename_beamform(array_str, time_start, freq_min, freq_max)

beamform_output = load_beamform_output(array_str, time_start, time_stop, 
                                                 freq_min=freq_min, freq_max=freq_max, 
                                                 slow_min=1/0.40, slow_max=1/0.30)


plot_backaz(beamform_output, settings.path_home, 
                        f"{array_str} Array, Filtered {freq_min} to {freq_max} Hz", 
                        file_str)


