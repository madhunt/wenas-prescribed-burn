#!/usr/bin/python3

import numpy as np
import pandas as pd
import itertools, datetime, pytz, os

import settings
import utils.load_data as load_data
settings.set_paths('laptop')

def main():
    time_start = datetime.datetime(2025, 10, 10, 16, 0, 0, tzinfo=pytz.timezone('UTC'))
    time_stop = datetime.datetime(2025, 10, 11, 5, 0, 0, tzinfo=pytz.timezone('UTC'))
    array_list = ['NW', 'NC', 'NE', 'SC']
    
    ints_all = calc_array_intersections(array_list, time_start, time_stop)

    source_loc = calc_mean_intersections(ints_all)

    format = '%Y-%m-%dT%H:%M'
    filename = f'source_locs_{time_start.strftime(format)}-{time_stop.strftime(format)}_4-8Hz' 
    source_loc.to_csv(os.path.join(settings.path_processed, "triang_source_locations", f"{filename}.csv"))
    source_loc.to_pickle(os.path.join(settings.path_processed, "triang_source_locations", f"{filename}.pkl"))
    return


def calc_array_intersections(array_list, time_start, time_stop):
    '''
    Docstring for calc_array_intersections
    
    :param array_list: Description
    :param time_start: Description
    :param time_stop: Description
    '''
    ints_all = pd.DataFrame()
    # loop through each pair of arrays
    for i, (arr1_str, arr2_str) in enumerate(itertools.combinations(array_list, 2)):
        # load in processed data and remove points with bad slownesses
        arr1 = load_data.load_beamform_output(arr1_str, 
                                              time_start=time_start, 
                                              time_stop=time_stop, 
                                              freq_min=4, freq_max=8, 
                                              slow_min=2.5, slow_max=3.5)
        arr2 = load_data.load_beamform_output(arr2_str, 
                                              time_start=time_start, 
                                              time_stop=time_stop, 
                                              freq_min=4, freq_max=8, 
                                              slow_min=2.5, slow_max=3.5)
        # find center of each array (easting, northing)
        p1 = array_coords_avg(array_str=arr1_str)
        p2 = array_coords_avg(array_str=arr2_str)

        # calculate intersections between this pair of arrays
        ints_pair = pd.DataFrame()

        ints_pair['result'] = arr1['Backaz'].combine(arr2['Backaz'], 
                                                     (lambda a1, a2: calc_intersection(p1, a1, p2, a2)))
        # format nicely in df
        ints_pair['Easting'] = [x[0] for x in ints_pair['result']]
        ints_pair['Northing'] = [x[1] for x in ints_pair['result']]
        ints_pair = ints_pair.drop('result', axis=1)      # clean up temp column

        # save intersection points in larger df
        ints_all[f'{arr1_str} {arr2_str} Easting'] = ints_pair['Easting']
        ints_all[f'{arr1_str} {arr2_str} Northing'] = ints_pair['Northing']

    return ints_all


def calc_mean_intersections(ints_all):
    '''
    Docstring for calc_mean_intersections
    
    :param ints_all: Description
    '''
    # calculate trimmed mean of intersections for each step in time
    source_loc = pd.DataFrame(index=ints_all.index,
                              columns=['Easting', 'Northing', 'Total Ints', 
                                       'NW Ints', 'NC Ints', 'NE Ints', 'SC Ints'])
    for t, row in ints_all.iterrows():
        #valid_ints = len(row[row.notna()])
        eastings = row.filter(like='Easting')#.to_numpy()
        northings = row.filter(like='Northing')#.to_numpy()

        # remove outliers
        outlier_cond = lambda ints: ((ints > (np.nanmean(ints) + 2*np.nanstd(ints))) | 
                                       (ints < (np.nanmean(ints) - 2*np.nanstd(ints))))
        inlier_east = ~outlier_cond(eastings).to_numpy()
        inlier_north = ~outlier_cond(northings).to_numpy()
        eastings = eastings[inlier_east & inlier_north]
        northings = northings[inlier_east & inlier_north]

        # calculate trimmed mean
        eastings_mean = np.nanmean(eastings)
        northings_mean = np.nanmean(northings)

        # calculate number of non-nan and non-outlier intersections
        num_ints = np.sum(~np.isnan(eastings))

        # store number of intersections with each array
        for array_str in ['NW', 'NC', 'NE', 'SC']:
            source_loc[f'{array_str} Ints'][t] = np.sum(~np.isnan(eastings.filter(like=array_str)))

        # store in dataframe
        source_loc['Easting'][t] = eastings_mean
        source_loc['Northing'][t] = northings_mean
        source_loc['Total Ints'][t] = num_ints
    return source_loc


def array_coords_avg(array_str, remove_bad=True):
    '''
    Docstring for array_coords_avg
    
    :param array_str: Description
    :param remove_bad: Description
    '''
    coords = pd.read_csv(settings.path_coords)
    coords = coords[coords['Network'] == array_str]

    if remove_bad:
        bad_stations = ['123', '005', '086', '368']
        idxs_bad = coords['SN'].str.contains('|'.join(bad_stations))
        coords = coords[~idxs_bad]

    # calculate mean easting/northing
    easting = coords['Easting'].mean()
    northing = coords['Northing'].mean()
    return easting, northing


def calc_intersection(p1, a1, p2, a2):
    '''
    Calculates intersection point between two rays, given starting points and azimuths (from N).
    INPUTS
        p1      : np array, 2x1 : Coordinates for start point of ray 1.
        a1      : float         : Azimuth of ray 1 direction in degrees (clockwise from North).
        p2      : np array, 2x1 : Coordinates for start point of ray 2.
        a2      : float         : Azimuth of ray 2 direction in degrees (clockwise from North).
    RETURNS
        int_pt : np array, 2x1 : Coordinates of intersection point (x, y). 
            Returns [NaN, NaN] if there is no intersection; e.g. if intersection 
            occurs "behind" ray start points or if rays are parallel. 
    '''
    if np.isnan(a1) or np.isnan(a2):
        return np.array([np.nan, np.nan])

    # create matrix of direction unit vectors
    D = np.array([[np.sin(np.radians(a1)), -np.sin(np.radians(a2))],
                      [np.cos(np.radians(a1)), -np.cos(np.radians(a2))]])

    # create vector of difference in start coords (p2x-p1x, p2y-p1y)
    P = np.array([p2[0] - p1[0],
                  p2[1] - p1[1]])

    # solve system of equations Dt=P
    try:
        t = np.linalg.solve(D, P)
    except:
        # matrix is singular (rays are parallel)
        return np.array([np.nan, np.nan])

    # see if intersection point is actually along rays
    if t[0]<0 or t[1]<0:
        # if intersection is "behind" rays, return nans
        return np.array([np.nan, np.nan])
        
    # calculate intersection point
    int_pt = np.array([p1[0]+D[0,0]*t[0],
                       p1[1]+D[1,0]*t[0]])
    return int_pt


if __name__ == "__main__":
    main()