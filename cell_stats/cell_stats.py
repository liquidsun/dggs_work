import sys
import os
import geopandas as gpd
from shapely.geometry import Polygon, shape
import argparse
import glob
import pandas as pd
import json
from multiprocessing import Pool
sys.path.append('..')
from h3_helper import *
from s2_helper import *
import numpy as np
import functools
from pandas.core.common import flatten
from dggrid4py import DGGRIDv7, Dggs, dgselect, dggs_types
from rhealpixdggs.dggs import *

# rhealpix functions
# TODO move to separate module

def get_rhpix_cells(res, extent=None):
    rdggs = WGS84_003
    if extent:
        nw = (extent[1], extent[2])
        se = (extent[3], extent[0])
        set_hex = list(flatten(rdggs.cells_from_region(res, nw, se, plane=False)))
    else:    
        set_hex = [x for x in rdggs.grid(res)]

    df = pd.DataFrame({"cell_id": set_hex})
    
    return df

def lonlat_to_latlon(lonlat_array):
    latlon_array = []
    for vertex in lonlat_array:
        latlon_array.append((vertex[1],vertex[0]))
    return latlon_array

def create_rhpix_geometry(df):
    gdf = gpd.GeoDataFrame(df.copy())
    gdf['geometry'] = gdf['cell_id'].apply(lambda x: Polygon(lonlat_to_latlon(x.boundary(n=2,plane=False))))
    gdf.crs = 'EPSG:4326'
    gdf['cell_id'] = gdf['cell_id'].apply(lambda x: str(x))
    return gdf

# rhealpix bloc finished 


def timer(func):

    @functools.wraps(func)
    def wrapper_timer(*args,**kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"finished{func.__name__!r} in {run_time:.4f} sec")
        return value
    return wrapper_timer


def get_conf(conf_path):
    with open(conf_path, 'r') as f:
        config = json.load(f)
    return config

@timer
def get_cells_area_stats(df, res):
    print(df.head())
    area_min = df['area'].min()
    area_max = df['area'].max()
    area_std = df['area'].std()
    area_mean = df['area'].mean()
    num_cells = len(df)
    
    stats_pd = pd.DataFrame({'resolution':[res],'min_area':[area_min],'max_area':[area_max],\
                             'std':[area_std],'mean':[area_mean],'num_cells':[num_cells]})
    
    return stats_pd

@timer
def get_cells_area(gdf,crs):
    
    gdf = gdf.to_crs(crs)
    gdf['area'] = gdf['geometry'].area
    
    return gdf

@timer
def create_cells(dggs, resolution, dggrid, extent=None):

    if dggs[0] == 'h3':
        if extent:
            df = get_h3_cells(resolution,extent['geometry'])
        else:
            df = get_h3_cells(resolution,extent)
        
        gdf = create_h3_geometry(df)

    elif dggs[0] == 's2':
        if extent:
            df = get_s2_cells(resolution,extent['bbox'])
        else:
            df = get_s2_cells(resolution,extent)
      
        gdf = create_s2_geometry(df)

    elif dggs[0] == 'DGGRID':
        if extent:
            extent = shape(extent['geometry'])
            gdf = dggrid.grid_cell_polygons_for_extent(dggs[1], resolution, clip_geom=extent)
            gdf.crs = 'EPSG:4326'
        else:
            gdf = dggrid.grid_cell_polygons_for_extent(dggs[1], resolution)
    
    elif dggs[0] == 'rhpix':
        if extent:
            df = get_rhpix_cells(resolution,extent['bbox'])
        else:
            df = get_rhpix_cells(resolution,extent)
      
        gdf = create_rhpix_geometry(df)

    return gdf


def cell_stats_parallel(func, geom_df, params, cores):

    extents = np.array_split([x for x in geom_df['geometry'].__geo_interface__['features']], cores)
    par_params = []
    for extent in extents:
        par = params.copy()
        par.append(extent)
        par_params.append(par)

    with Pool(cores) as p:
        ret_list = p.map(func, par_params)

    ret_list = pd.concat(ret_list)
    
    ret_list = ret_list.groupby('resolution').agg(resolution=('resolution','first'), min_area=('min_area','min'), max_area=('max_area','max'),
                                                  std=('std','mean'),mean =('mean','mean'),num_cells=('num_cells','sum'))

    return ret_list


def create_cell_stats_df(params):

    stats = []
    if len(params) == 4:
        gdf = create_cells(params[0], params[1], params[3])
        gdf = get_cells_area(gdf, params[2])
        stats.append(get_cells_area_stats(gdf,params[1]))
    
    else:
        for extent in params[4]:
            gdf = create_cells(params[0], params[1], params[3], extent)
            gdf = get_cells_area(gdf, params[2])
            stats.append(get_cells_area_stats(gdf,params[1]))

    return pd.concat(stats)



def main():

    #Getting the directory of the script
    script_dir = os.getcwd()

    #Parsing arguments
    parser = argparse.ArgumentParser()

    #Specifies the config file
    parser.add_argument('-config', action='store', dest='config', default=os.path.join(script_dir,'config.json'),
                        help='Specifies the folder where config files are')

    #Specify base output name
    parser.add_argument('-name', action='store', dest='result_name', default='cell_stats',
                        help='Name for export to database or csv')

    #Specify number of paralell processes
    parser.add_argument('-cpu', action='store', dest='cpu', default='3',
                        help='Number of processes to run')

    #Specify output folder
    parser.add_argument('-out', action='store', dest='out', default=os.path.join(script_dir,'results'),
                        help='Number of processes to run')

    #Specify sample polygons path
    parser.add_argument('-sample', action='store', dest='sample', default=os.path.join(script_dir,'sample_polygons.geojson'),
                        help='Number of processes to run')

    #Specify sample polygons path
    parser.add_argument('-dggrid', action='store', dest='dggrid', default=os.path.join(script_dir,'dggrid'),
                        help='Number of processes to run')

    #Get the results of argument parsing
    results = parser.parse_args()

    #Store the results of argument parsing

    results_path = results.out
    config_path = results.config
    default_name = results.result_name
    cpus = int(results.cpu)
    sample = results.sample
    dggrid_exec = results.dggrid

    sample_polygons = gpd.read_file(sample,driver='GeoJSON')
    config = get_conf(config_path)

    dggrid_work_dir = os.path.join(script_dir,'dggrid_workdir')

    if not os.path.exists(dggrid_work_dir):
        os.makedirs(dggrid_work_dir)

    for dggs in config['dggss']:
        print(f"Start processing {dggs['name']}")
        if dggs['name'][0] == 'DGGRID':
            dggrid_instance = DGGRIDv7(executable=dggrid_exec, working_dir=dggrid_work_dir, capture_logs=False, silent=False)
        else:
            dggrid_instance = None
        stats_df_global = []
        for res in dggs['global_res']:
            print(f"Start processing global resolution {res}")
            stats_df_global.append(create_cell_stats_df([dggs['name'], res, dggs['proj'], dggrid_instance]))

        for res in dggs['sample_res']:
            print(f"Start processing local resolution {res}")
            stats_df_sample = cell_stats_parallel(create_cell_stats_df, sample_polygons,[dggs['name'], res, dggs['proj'], dggrid_instance], cpus)
        
        final_stats = pd.concat([pd.concat(stats_df_global),stats_df_sample])

        if not os.path.exists(results_path):
            os.makedirs(results_path)
        if len(dggs['name']) == 1:
            name = dggs['name'][0]
        else:
            name = dggs['name'][0] + '_' + dggs['name'][1]

        res_file_name = os.path.join(results_path, default_name + f"_{name}")

        final_stats.to_csv(res_file_name, index=False)

if __name__ == "__main__":
    main()