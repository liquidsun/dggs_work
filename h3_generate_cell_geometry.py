import h3
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import transform
from geoalchemy2 import WKTElement, Geometry
import time
import pyproj
from functools import partial


def create_h3_geom_cells(extent_json, resolutions, table,  db_engine, export_type):
    """Create geometry for h3 cells in given extent for given resolutions levels

    Parameters:
    extentJSON (geojson): extent for covering with tiles
    resolutions([int]): array of integer h3 resolution levels
    table(string): table name for postgres database
    db_engine (sqlalchemy.engine): sqlalchemy database engine
    export_type(string): where to export 'geojson' or 'postgres'

    Returns:
    none
   """
    for res in resolutions:

        print('start caclulating resolution {} {}'.format(res, time.asctime(time.localtime(time.time()))))
        set_hex = list(h3.polyfill(extent_json['features'][0]["geometry"], res=res))

        print('finish caclulating resolution {} {}'.format(res, time.asctime(time.localtime(time.time()))))

        if export_type == 'postgres':
            df = pd.DataFrame({"hex_id": set_hex})
            df['geom'] = df.hex_id.apply(lambda x:
                                         WKTElement(Polygon(h3.h3_to_geo_boundary(x)).wkt, srid=4326))
            print('finish caclulating geometry {} {}'.format(res, time.asctime(time.localtime(time.time()))))
            df.to_sql(table + str(res), db_engine, if_exists='replace', index=True, schema='public',
                      dtype={'geom': Geometry('POLYGON', srid=4326)})
            print('finish import to db {} {}'.format(res, time.asctime(time.localtime(time.time()))))

        elif export_type == 'geojson':
            proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'),
                           pyproj.Proj('proj=sinu'))
            gdf = gpd.GeoDataFrame({"hex_id": set_hex})

            gdf['geometry'] = gdf.hex_id.apply(lambda x: Polygon(h3.h3_to_geo_boundary(x, True)))
            gdf['area'] = gdf.geometry.apply(lambda x: transform(proj, x).area)
            gdf.to_file("{}{}.geojson".format(table, res), driver='GeoJSON')
            print('finish import to geojson {} {}'.format(res, time.asctime(time.localtime(time.time()))))


def create_h3_geom_cells_global(resolutions, table, export_type, db_engine=''):
    """Create geometry for h3 cells globally for given resolutions

        Parameters:
        db_engine (sqlalchemy.engine): sqlalchemy database engine
        resolutions([int]): array of integer h3 resolution levels
        table(string): table name for postgres database
        export_type(string): where to export 'geojson' or 'postgres'

        Returns:
        none
    """
    for res in resolutions:
        set_hex_0 = list(h3.h3_get_res_0_indexes())
        set_hex = []
        if res == 0:
            set_hex = set_hex_0
        else:
            for i in set_hex_0:
                set_hex.extend(list(h3.h3_to_children(i, res)))
        if export_type == 'postgres':
            df = pd.DataFrame({"hex_id": set_hex})
            df['geom'] = df.hex_id.apply(lambda x:
                                         WKTElement(Polygon(h3.h3_to_geo_boundary(x)).wkt, srid=4326))
            print('finish caclulating geometry {} {}'.format(res, time.asctime(time.localtime(time.time()))))
            df.to_sql(table + str(res), db_engine, if_exists='replace', index=True, schema='public',
                      dtype={'geom': Geometry('POLYGON', srid=4326)})
            print('finish import to db {} {}'.format(res, time.asctime(time.localtime(time.time()))))

        elif export_type == 'geojson':
            proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'),
                           pyproj.Proj('proj=isea'))
            gdf = gpd.GeoDataFrame({"hex_id": set_hex})
            gdf['geometry'] = gdf.hex_id.apply(lambda x: Polygon(h3.h3_to_geo_boundary(x, True)))
            print('finish caclulating geometry {} {}'.format(res, time.asctime(time.localtime(time.time()))))
            gdf['area'] = gdf.geometry.apply(lambda x: transform(proj, x).area)
            gdf.to_file("{}{}.geojson".format(table, res), driver='GeoJSON')
            print('finish import to geojson {} {}'.format(res, time.asctime(time.localtime(time.time()))))





