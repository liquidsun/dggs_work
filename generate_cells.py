from sqlalchemy.engine.url import URL
from sqlalchemy import create_engine
from h3_helper import *
from s2_generate_cell_geometry import *
from EAGGR_generate_cell_geometry import *
from shapely.geometry import box



HOST = 'localhost'
DB = 'dggs'
USER = 'Iam'
PORT = 5432
PWD = 'postgres'

workDataDir = r"D:\UTCloud\DGGS\working\Python\workingdata"
# Database info
db_url = URL(drivername='postgresql+psycopg2', host=HOST, database=DB,
             username=USER, port=PORT, password=PWD)

# Create engine
engine = create_engine(db_url)

# read extent file
gdf = gpd.read_file(workDataDir + "working.gpkg", layer='dev_extent2')

# convert extent to geoJSON
extentJSON = json.loads((gdf.to_json()))

# create geojson extent
extentJSON = gpd.GeoSeries(box(12.2,47.3,13.7,48.2)).__geo_interface__

# create local h3 geometry and load it to database
create_h3_geom_cells(extentJSON, [9], 'h3_cells', engine, 'postgres')

# create global h3 geometry to geojson
create_h3_geom_cells_global([x for x in range(5)], 'h3_level', 'geojson')

# create global s2 geometry
create_s2_geom_cells(box(-180, -80, 180, 80), [x for x in range(8)])

# create global eaggrt geometry
create_eaggrt_cells_global([x for x in range(7)], 'data/eaggrt_level')

# create global eaggrh geometry
create_eaggrh_cells_global([x for x in range(5)], 'data/eaggrh_level')

