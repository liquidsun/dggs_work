

from io import StringIO
from h3 import h3
import pandas as pd
import geopandas as gpd
from sqlalchemy.engine.url import URL
from sqlalchemy import create_engine
import gdal
import csv
# import matplotlib
import time
import geojson


# Functions

def df_to_sql(db_engine, df, table, if_exists='replace', sep='\t', encoding='utf8'):
    # Create Table
    df[:0].to_sql(table, db_engine, if_exists=if_exists, index=False)
    print("{} table prepared".format(table))
    # Prepare data
    output = StringIO()
    df.to_csv(output, sep=sep, index=False, header=False, escapechar="\\", encoding=encoding, quoting=csv.QUOTE_NONE)
    output.seek(0)
    # Insert data
    connection = db_engine.raw_connection()
    print("Opened database successfully")
    insert_start_time = time.time()
    cursor = connection.cursor()
    cursor.copy_from(output, table, sep=sep, null='')
    connection.commit()
    print("insert finish %s seconds ---" % (time.time() - insert_start_time))
    connection.close()


def raster_centroids_to_h3(work_dir, raster, hex_scale):
    hex_col = 'hex' + str(hex_scale)
    # convert .tif to .xyz
    gdal.Translate(work_dir + "dev_elevationXYZ.xyz", raster, format="XYZ")

    # read .xvz to data frame
    raster_reading_time = time.time()
    df = pd.read_csv(work_dir + 'dev_elevationXYZ.xyz', sep=' ', header=None, names=['lng', 'lat', 'elevation'])
    print("Read raster done {}".format(time.time() - raster_reading_time))

    # Get h3 indexes from coordinate pairs
    geo_to_hex_time = time.time()
    df[hex_col] = df.apply(lambda x: h3.geo_to_h3(x.lng, x.lat, hex_scale), 1)
    print("Geo_to_hex_time{}".format(time.time() - geo_to_hex_time))

    # group equal indexes with calculation of mean value
    df_dem = df.groupby(hex_col)['elevation'].mean().to_frame('elevation').reset_index()
    return df_dem


def fill_h3_with_raster_val(raster, hex_scale, extent_JSON, name):
    set_hexagons = list(h3.polyfill(extent_JSON['features'][0]["geometry"], res=hex_scale))
    rs = gdal.Open(raster)
    gt = rs.GetGeoTransform()
    rb = rs.GetRasterBand(1)
    df = pd.DataFrame({"hex_id": set_hexagons})
    df[name] = df.hex_id.apply(lambda x: rb.ReadAsArray(int((h3.h3_to_geo(x)[0] - gt[0]) / gt[1]),
                                                        int((h3.h3_to_geo(x)[1] - gt[3]) / gt[5]), 1, 1)[0][0])
    return df


def hex_h3_downsampling(df, hex_col, metric_col, coarse_hex_scale):
    df_coarse = df.copy()
    coarse_hex_col = 'hex{}'.format(coarse_hex_scale)
    df_coarse[coarse_hex_col] = df_coarse[hex_col].apply(lambda x: h3.h3_to_parent(x, coarse_hex_scale))
    dfc = df_coarse.groupby([coarse_hex_col])[[metric_col, ]].mean().reset_index()
    return dfc


# Create db engine
HOST = 'venuz'
DB = 'DGGS'
USER = 'postgres'
PORT = 5432
PWD = 'postgres'
db_url = URL(drivername='postgresql+psycopg2', host=HOST, database=DB,
             username=USER, port=PORT, password=PWD)
engine = create_engine(db_url)

workDataDir = r"./data/workingdata/"
hex_size = 11

#  Read extent for h3 cells
gdf = gpd.read_file(workDataDir + "working.gpkg", layer='dev_extent2')
extentJSON = geojson.loads((gdf.to_json()))

#  Read raster to load
raster = r"./data/EestiClipped/small_extent/MeritDEM90m.tif"

#  Load raster values to h3 cells, create geodataframe
elev_df = fill_h3_with_raster_val(raster, hex_size, extentJSON, 'elevation')

#  Load geodataframe to Postgres
df_to_sql(engine, elev_df, 'h3_' + str(hex_size) + '_elev_dev')
