from eaggr.eaggr import Eaggr
from eaggr.enums.shape_string_format import ShapeStringFormat
from eaggr.shapes.dggs_cell import DggsCell
import geopandas as gpd
import pandas as pd
from eaggr.enums.model import Model
import json
from shapely.geometry import Polygon
import pyproj
from functools import partial
from shapely.ops import transform


def get_eaggr_indexes_at_level(df_level_0, resolution, dggs):
    if resolution == 0:
        return df_level_0
    else:
        next_res = gpd.GeoDataFrame(df_level_0['cell'].apply(lambda x: pd.Series(dggs.get_dggs_cell_children(x))) \
                                    .stack().to_frame('cell').reset_index(1, drop=True).reset_index(drop=True))
        return get_eaggr_indexes_at_level(next_res, resolution - 1, dggs)


def create_eaggrh_cells_global(resolutions, table):
    dggs = Eaggr(Model.ISEA3H)
    cell_ids_0 = ['03000,0', '13000,0', '14000,0', '09000,0', '17000,0', '15000,0', '10000,0', '01000,0', '12000,0',
                  '16000,0',
                  '07000,0', '08000,0', '11000,0', '19000,0', '04000,0', '05000,0', '06000,0', '00000,0', '18000,0',
                  '02000,0']
    gdf_level_0 = gpd.GeoDataFrame()
    gdf_level_0['cell'] = pd.Series(list(map(lambda x: DggsCell(x), cell_ids_0)))
    proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'),
                   pyproj.Proj('proj=isea'))

    for r in resolutions:
        gdf_level_ = get_eaggr_indexes_at_level(gdf_level_0, r, dggs).copy()
        gdf_level_['cell_id'] = gdf_level_['cell'].apply(lambda x: x.get_cell_id())
        gdf_level_.drop_duplicates('cell_id')
        gdf_level_['geojson'] = gdf_level_['cell'].apply(lambda x: json.loads
        (dggs.convert_dggs_cell_outline_to_shape_string(x, ShapeStringFormat.GEO_JSON)))

        if r == 0:
            gdf_level_['geometry'] = gdf_level_['geojson'].apply(
                lambda x: Polygon([x['coordinates'][0][0], x['coordinates'][0][1],
                                   x['coordinates'][0][2]]))
        else:
            gdf_level_['geometry'] = gdf_level_['geojson'].apply(
                lambda x: Polygon([x['coordinates'][0][0], x['coordinates'][0][1],
                                   x['coordinates'][0][2], x['coordinates'][0][3], x['coordinates'][0][4],
                                   x['coordinates'][0][5]]))

        gdf_level_.drop(columns=['cell', 'geojson'], inplace=True, index=1)

        gdf_level_['area'] = gdf_level_['geometry'].apply(lambda x: transform(proj, x).area)

        gdf_level_.to_file("{}{}.geojson".format(table, r), driver='GeoJSON')


def create_eaggrt_cells_global(resolutions, table):
    dggs = Eaggr(Model.ISEA4T)
    cell_ids_0 = list(map(lambda x: '0' + x, list(map(str, [x for x in range(10)]))))
    cell_ids_0.extend(list(map(str, [x for x in range(10, 20, 1)])))
    gdf_level_0 = gpd.GeoDataFrame()
    gdf_level_0['cell'] = pd.Series(list(map(lambda x: DggsCell(x), cell_ids_0)))

    proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'),
                   pyproj.Proj('proj=isea'))

    for r in resolutions:
        gdf_level_ = get_eaggr_indexes_at_level(gdf_level_0, r, dggs).copy()
        gdf_level_['cell_id'] = gdf_level_['cell'].apply(lambda x: x.get_cell_id())

        gdf_level_['geojson'] = gdf_level_['cell'].apply(lambda x: json.loads
        (dggs.convert_dggs_cell_outline_to_shape_string(x, ShapeStringFormat.GEO_JSON)))

        gdf_level_['geometry'] = gdf_level_['geojson'].apply(
            lambda x: Polygon([x['coordinates'][0][0], x['coordinates'][0][1],
                               x['coordinates'][0][2]]))

        gdf_level_.drop(columns=['cell', 'geojson'], inplace=True, index=1)

        gdf_level_['area'] = gdf_level_['geometry'].apply(lambda x: transform(proj, x).area)

        gdf_level_.to_file("{}{}.geojson".format(table, r), driver='GeoJSON')


create_eaggrh_cells_global([x for x in range(5)], 'data/eaggrh_level')
