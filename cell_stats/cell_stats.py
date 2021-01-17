import sys
import os
import geopandas as gpd
from shapely.geometry import Polygon, shape
from shapely.ops import transform
import pyproj
import argparse
import glob
import pandas as pd
import json
import math
from multiprocessing import Pool
sys.path.append('..')
from h3_helper import *
from s2_helper import *
from rhealpix_helper import *
import numpy as np
import functools
from dggrid4py import DGGRIDv7, Dggs, dgselect, dggs_types


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


def get_area_perimeter_from_lambert(geom):
    '''Area from cell's lambert azimutal projection'''
    try:
        if (-180 <= geom.centroid.x <= 180) and (-90 <= geom.centroid.y <= 90):
            proj_str = f"+proj=laea +lat_0={geom.centroid.y} +lon_0={geom.centroid.x}"
            project = pyproj.Transformer.from_crs('EPSG:4236', proj_str, always_xy=True).transform
            perimeter = transform(project, geom).length
            area = transform(project, geom).area
        else:
            print(f'invalid centroid {geom.centroid}')
            perimeter = None
            area = None
    except:
        print(f'invalid centroid None')
        perimeter = None
        area = None
    return area, perimeter

def zsc_calculation(df):
    zsc = math.sqrt(4*math.pi*df['area'] - math.pow(df['area'],2)/math.pow(6378137,2))/df['perimeter']
    return zsc

def check_crossing(lon1: float, lon2: float, validate: bool = True):
    """
    Assuming a minimum travel distance between two provided longitude coordinates,
    checks if the 180th meridian (antimeridian) is crossed.
    """
    if validate and any(abs(x) > 180.0 for x in [lon1, lon2]):
        raise ValueError("longitudes must be in degrees [-180.0, 180.0]")   
    return abs(lon2 - lon1) > 180.0

def check_for_geom(geom):
    crossed = False
    p_init = geom.exterior.coords[0]

    for p in range(1, len(geom.exterior.coords)):
        px = geom.exterior.coords[p]
        # print(px)

        if check_crossing(p_init[0], px[0]):
            crossed = True
        p_init = px
    
    return crossed

@timer
def get_cells_area_stats(df, res):
    
    # Filter out invalid geometry
    try:
        df['crossed'] = df['geometry'].apply(check_for_geom)
        date_line_cross_error_cells = len(df[df['crossed']])
        df = df[~df['crossed']]
        other_geom_anomalies = len(df[(df['area']<df['area'].quantile(0.005))&(df['area']>df['area'].quantile(0.995))])
        df = df[(df['area']>df['area'].quantile(0.005))&(df['area']<df['area'].quantile(0.995))]
        df['std_area'] = df['area']/df['area'].mean() 
        df['zsc'] = df.apply(zsc_calculation,axis=1)
        area_min = df['area'].min()
        area_max = df['area'].max()
        area_std = df['area'].std()
        area_mean = df['area'].mean()
        std_area_std = df['std_area'].std()
        std_area_range = df['std_area'].max() - df['std_area'].min()
        zsc_std = df['zsc'].std()
        zsc_std_range = df['zsc'].max() - df['zsc'].min()
        num_cells = len(df)
        
        stats_pd = pd.DataFrame({'resolution':[res],'min_area':[area_min],'max_area':[area_max],\
                                'std':[area_std],'mean':[area_mean],'num_cells':[num_cells], 'std_area_std':[std_area_std],\
                                'std_area_range':[std_area_range], 'zsc_std':[zsc_std], 'zsc_std_range':[zsc_std_range],\
                                'date_line_cross_error_cells':[date_line_cross_error_cells],'other_geom_anomalies':other_geom_anomalies})
    
        return stats_pd
    except: 
        print(len(df))
        return None


def get_lambert_area(geom):
    '''Area from cell's lambert azimutal projection'''

    if (-180 <= geom.centroid.x <= 180) and (-90 <= geom.centroid.y <= 90):
        proj_str = f"+proj=laea +lat_0={geom.centroid.y} +lon_0={geom.centroid.x}"
        project = pyproj.Transformer.from_crs('EPSG:4236', proj_str, always_xy=True).transform
        area = transform(project, geom).area
    else:
        print(f'invalid centroid {geom.centroid}')
        area = None
    return area


def convert_wgs_to_utm(lon, lat):
    '''Get UTM zone from lat lon''' 

    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
    else:
        epsg_code = '327' + utm_band
    return epsg_code


def get_utm_area(geom):
    '''Area from cell's UTM zone'''
    if (-180 <= geom.centroid.x <= 180) and (-90 <= geom.centroid.y <= 90):
        utm_epsg = convert_wgs_to_utm(geom.centroid.x, geom.centroid.y) 
        project = pyproj.Transformer.from_crs('EPSG:4236', f'EPSG:{utm_epsg}', always_xy=True).transform
        area = transform(project, geom).area
    else:
        print(f'invalid centroid {geom.centroid}')
        area = None
    return area


@timer
def get_cells_area(gdf,crs):
    '''Get cells area for crs'''

    if crs =='UTM':
        gdf['area'] = gdf['geometry'].apply(get_utm_area)

    elif crs =='LAEA':
        gdf['area'],gdf['perimeter'] = zip(*gdf['geometry'].apply(get_area_perimeter_from_lambert))
        
    else:
        gdf = gdf.to_crs(crs)
        gdf['area'] = gdf['geometry'].area
    return gdf

@timer
def create_cells(dggs, resolution, dggrid, extent=None):
    '''Creates cells for given DGGS on given resolution'''

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
    ''' Split sample areas stats calculations between cores '''
    extents = np.array_split([x for x in geom_df['geometry'].__geo_interface__['features']], cores)
    par_params = []
    for extent in extents:
        par = params.copy()
        par.append(extent)
        par_params.append(par)

    with Pool(cores) as p:
        ret_list = p.map(func, par_params)

    ret_list = pd.concat(ret_list)
    
    ret_list = ret_list.groupby('resolution').agg(resolution=('resolution','first'), min_area=('min_area','min'),std_area_std=('std_area_std','mean'), max_area=('max_area','max'),
                                                  std=('std','mean'),mean = ('mean','mean'),num_cells=('num_cells','sum'),std_area_range=('std_area_range','mean'),
                                                  zsc_std=('zsc_std','mean'),zsc_std_range = ('zsc_std_range','mean'),date_line_cross_error_cells=('date_line_cross_error_cells','sum'),
                                                  other_geom_anomalies=('other_geom_anomalies','sum')
                                                  )

    return ret_list



def create_cell_stats_df(params):
    '''Create cell stats. Global coverage or extnet coverage'''
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

        stats_df_sample = []
        for res in dggs['sample_res']:
            print(f"Start processing local resolution {res}")
            stats_df_sample.append(cell_stats_parallel(create_cell_stats_df, sample_polygons,[dggs['name'], res, dggs['proj'], dggrid_instance], cpus))
        
        if len(stats_df_sample) > 0:
            final_stats = pd.concat([pd.concat(stats_df_global),pd.concat(stats_df_sample)])
        else:
            final_stats = pd.concat([pd.concat(stats_df_global)])

        
        if len(dggs['name'])>1:
            final_stats['dggs'] = dggs['name'][1]
        else:
            final_stats['dggs'] = dggs['name'][0]
            
        final_stats['proj'] = dggs['proj']
        if not os.path.exists(results_path):
            os.makedirs(results_path)
        if len(dggs['name']) == 1:
            name = dggs['name'][0]
        else:
            name = dggs['name'][0] + '_' + dggs['name'][1]

        res_file_name = os.path.join(results_path, default_name + f"_{name}.csv")

        final_stats.to_csv(res_file_name, index=False)

if __name__ == "__main__":
    main()