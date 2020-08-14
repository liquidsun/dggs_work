import pandas as pd
import geopandas as gpd
import time
import geojson
from pyproj import Transformer
import rasterio
from eaggr.eaggr import Eaggr
from eaggr.shapes.lat_long_point import LatLongPoint
from eaggr.enums.shape_string_format import ShapeStringFormat
from eaggr.shapes.dggs_cell import DggsCell
import eaggr.shapes
from eaggr.enums.model import Model
from shapely.geometry import box, Polygon, Point
import numpy as np
import json
from shapely.ops import transform


def generate_grid(extent, grid_size):
    """Generate sampling grid points in given extent with given number of points

    Parameters:
    extent (shapely box): extent to cover with sampling points
    grid_size (integer): number of sampling points

    Returns:
    Pandas dataframe
   """

    transformer = Transformer.from_crs("epsg:4326", 'epsg:3301')

    # convert to planar
    exten_plain = box(transformer.transform(extent.bounds[1], extent.bounds[0])[1],
                      transformer.transform(extent.bounds[1], extent.bounds[0])[0],
                      transformer.transform(extent.bounds[3], extent.bounds[2])[1],
                      transformer.transform(extent.bounds[3], extent.bounds[2])[0])

    # calculate parameters for the grid
    x_min = exten_plain.bounds[0]
    y_min = exten_plain.bounds[1]
    x_size = exten_plain.bounds[2] - exten_plain.bounds[0]
    y_size = exten_plain.bounds[3] - exten_plain.bounds[1]

    grid_x_res = x_size / grid_size
    grid_y_res = y_size / grid_size

    # calculate x,y of grid with real coordinates
    x_grid = [(x * grid_x_res) + x_min for x in range(grid_size)]
    y_grid = [(y * grid_y_res) + y_min for y in range(grid_size)]

    # calculate all x,y values of grid
    x_grid_coord, y_grid_coord = np.meshgrid(x_grid, y_grid)

    # create geoseries with grid points and assign projection
    int_grid = gpd.GeoSeries(map(Point, zip(x_grid_coord.flatten(), y_grid_coord.flatten())))
    int_grid.crs = "EPSG:3301"

    # transform back to lat/lon
    int_grid = int_grid.to_crs("EPSG:4326")
    int_grid = gpd.GeoDataFrame({'geometry': int_grid})
    return int_grid


def get_eaggr_indexes_at_level(df_level_0, resolution, dggs):
    if resolution == 0:
        return df_level_0
    else:
        next_res = gpd.GeoDataFrame(df_level_0['cell'].apply(lambda x: pd.Series(dggs.get_dggs_cell_children(x))) \
                                    .stack().to_frame('cell').reset_index(1, drop=True).reset_index(drop=True))
        return get_eaggr_indexes_at_level(next_res, resolution - 1, dggs)


def create_eaggrt_cells_global(resolutions, table):
    dggs = Eaggr(Model.ISEA4T)
    cell_ids_0 = list(map(lambda x: '0' + x, list(map(str, [x for x in range(10)]))))
    cell_ids_0.extend(list(map(str, [x for x in range(10, 20, 1)])))
    gdf_level_0 = gpd.GeoDataFrame()
    gdf_level_0['cell'] = pd.Series(list(map(lambda x: DggsCell(x), cell_ids_0)))

    transformer = Transformer.from_crs("epsg:4326", 'proj=isea')

    for r in resolutions:
        gdf_level_ = get_eaggr_indexes_at_level(gdf_level_0, r, dggs).copy()
        gdf_level_['cell_id'] = gdf_level_['cell'].apply(lambda x: x.get_cell_id())

        gdf_level_['geojson'] = gdf_level_['cell'].apply(lambda x: json.loads
        (dggs.convert_dggs_cell_outline_to_shape_string(x, ShapeStringFormat.GEO_JSON)))

        gdf_level_['geometry'] = gdf_level_['geojson'].apply(
            lambda x: Polygon([x['coordinates'][0][0], x['coordinates'][0][1],
                               x['coordinates'][0][2]]))

        gdf_level_.drop(columns=['cell', 'geojson'], inplace=True, index=1)

        gdf_level_['area'] = gdf_level_['geometry'].apply(lambda x: transform(transformer.transform, x).area)

        gdf_level_.to_file("{}{}.geojson".format(table, r), driver='GeoJSON')


def raster_to_eaggr_t(raster_path, value_name, sampling_grid_size, grid_point_accuracy, extent=None):
    """Load raster values into s2 dggs cells

    Parameters:
    raster (string): path to raster file for uploading
    value_name (string): name of a value to be uploaded
    sampling_grid_size (integer): number of sampling points
    grid_point_accuracy (integer): accuracy of sampling points in m
    extent (list): Extent as array of 2 lat lon pairs to get raster values for

    Returns:
    Pandas dataframe
   """

    dggs = Eaggr(Model.ISEA4T)
    # Open raster
    rs = rasterio.open(raster_path)

    # Get extent to fill with eaggr triangles
    if extent:
        extent = box(26.521263, 58.283163, 26.900083, 58.505180)
    else:
        extent = box(rs.bounds.bottom, rs.bounds.left, rs.bounds.top, rs.bounds.right)
    # Get resolution dict

    # Get two neighbour pixels in raster

    # Get pixel size from projected src

    # Get sampling grid

    eaggr_gdf = generate_grid(extent, sampling_grid_size)

    # Get raster band as np array
    raster_band_array = rs.read(1)

    # Get dggs resolution for raster pixel size

    # Fill dataframe with cell_ids  from cover with given resolution
    print(f"Start filling extent with eaggr indexes at resolution")
    eaggr_gdf['cell_id'] = eaggr_gdf['geometry'].apply(
        lambda x: dggs.convert_point_to_dggs_cell(LatLongPoint(latitude=x.y,
                                                               longitude=x.x,
                                                               accuracy=grid_point_accuracy)).get_cell_id())

    # Get raster values for each hex_id
    print(f"Start getting raster values for eaggr cells at resolution")
    eaggr_gdf[value_name] = eaggr_gdf['cell_id'].apply(
        lambda x: raster_band_array[rs.index(dggs.convert_dggs_cell_to_point(DggsCell(x)).get_longitude(),
                                             dggs.convert_dggs_cell_to_point(DggsCell(x)).get_latitude())])

    # Drop unnecessary columns
    eaggr_gdf = eaggr_gdf[eaggr_gdf[value_name] != rs.nodata]
    eaggr_gdf = eaggr_gdf.groupby('cell_id').first().reset_index()
    eaggr_gdf.drop('geometry', axis=1, inplace=True)
    return eaggr_gdf


def vector_to_eaggr_t(vector_path, value_name, sampling_grid_size, grid_point_accuracy, extent=None, layer=None):
    """Load vector values into s2 dggs cells

    Parameters:
    vector_path (string): path to vector file for uploading
    value_name (string): name of a vector attribute to be uploaded
    resolution (integer): s2 resolution to load vector values into
    extent (list): Extent as array of 2 lat lon pairs to get vector values for
    Returns:
    Pandas dataframe
   """

    dggs = Eaggr(Model.ISEA4T)

    # Open vector to geodataframe
    gdf = gpd.read_file(vector_path, layer)

    # Get extent to fill with s2 squares
    if extent == None:
        extent = box(gdf['geometry'].total_bounds[1], gdf['geometry'].total_bounds[0],
                     gdf['geometry'].total_bounds[2], gdf['geometry'].total_bounds[3])

    # Get sampling grid
    sampling_grid = generate_grid(extent, sampling_grid_size)

    # Fill dataframe with cell_ids  from cover with given resolution
    print(f"Start filling extent with eaggr indexes at resolution")
    sampling_grid['cell_id'] = sampling_grid['geometry'].apply(
        lambda x: dggs.convert_point_to_dggs_cell(LatLongPoint(latitude=x.y,
                                                               longitude=x.x,
                                                               accuracy=grid_point_accuracy)).get_cell_id())

    # Get hex centroids for points
    sampling_grid['geometry'] = sampling_grid['cell_id'].apply(
        lambda x: Point(dggs.convert_dggs_cell_to_point(DggsCell(x)).get_longitude(),
                        dggs.convert_dggs_cell_to_point(DggsCell(x)).get_latitude()))
    sampling_grid.set_crs('epsg:4326')

    # Spatial join hex centroids with gdf
    vector_eaggr = gpd.sjoin(sampling_grid, gdf)

    # Drop unnecessary fields
    vector_eaggr = vector_eaggr.groupby('cell_id').first().reset_index()
    vector_eaggr = vector_eaggr[['cell_id', value_name]]
    # vector_eaggr.drop('geometry',axis=1, inplace=True)

    return vector_eaggr


def cell_eaggr_t_downsampling(df, cell_id_col, metric_col, coarse_resolution, metric_type):
    """Aggregates a given attribute in eaggr cell to a given coarser resolution level

    Parameters:
    df (pandas dataframe): dataframe with s2 ids and attributes for aggregation
    cell_id_col (string): name of s2 id column
    metric_col (string): name of a column for aggreagation
    coarse_resolution (integer): Coarser s2 resoluiton for aggregation
    metric_type (string): attribute type (numerical, categorical)
    Returns:
    Pandas dataframe
   """

    dggs = Eaggr(Model.ISEA4T)
    df_coarse = df.copy()
    coarse_id_col = 'cell_id_{}'.format(coarse_resolution)
    df_coarse[coarse_id_col] = df_coarse[cell_id_col].apply(
        lambda x: dggs.get_dggs_cell_parents(DggsCell(x))[0].get_cell_id())

    if metric_type == 'numeric':
        dfc = df_coarse.groupby(coarse_id_col)[[metric_col]].mean().reset_index()
    elif metric_type == 'categorical':
        dfc = df_coarse.groupby([coarse_id_col, metric_col]).agg(count=(metric_col, 'count')).reset_index().sort_values(
            by=[coarse_id_col, metric_col, 'count']).groupby(coarse_id_col, as_index=False, sort=False).first()
        dfc.drop('count', axis=1, inplace=True)

    dfc.columns = [cell_id_col, metric_col]

    return dfc


def create_eaggr_geometry(df):

    dggs = Eaggr(Model.ISEA4T)
    df['geojson'] = df['cell_id'].apply(lambda x: json.loads
    (dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(x), ShapeStringFormat.GEO_JSON)))

    df['geometry'] = df['geojson'].apply(
        lambda x: Polygon([x['coordinates'][0][0], x['coordinates'][0][1],
                           x['coordinates'][0][2]]))

    df = gpd.GeoDataFrame(df)
    df.crs = "EPSG:4326"
    return df