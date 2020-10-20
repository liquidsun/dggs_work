import pywraps2 as s2
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point, box
from shapely.ops import transform
from pyproj import Transformer
import rasterio
from math import radians, sin, cos, asin, sqrt


def __haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r * 1000


def s2_geometry_from_cellid(cell_id):
    new_cell = s2.S2Cell(s2.S2CellId(cell_id))
    vertices = []
    for i in range(0, 4):
        vertex = new_cell.GetS2LatLngVertex(i)
        vertices.append((vertex.lng().degrees(),
                         vertex.lat().degrees()))
    geom = Polygon(vertices)
    return geom

def create_s2_geometry(df):
    gdf = gpd.GeoDataFrame(df)
    gdf['geometry'] = gdf['cell_id'].apply(lambda x: s2_geometry_from_cellid(x))
    gdf.crs = 'EPSG:4326'
    return gdf


def get_s2_cells(res, extent=None):
    """Get s2 cells for given resolution

    Parameters:
    res (int): S2 resolution 
    extent (list): Extent as array of 2 lon lat pairs to get raster values for
    Returns:
    Pandas dataframe
   """

    coverer = s2.S2RegionCoverer()
    if extent:
        coverer.set_fixed_level(res)
        region_rect = s2.S2LatLngRect(
        s2.S2LatLng_FromDegrees(extent[1], extent[0]),
        s2.S2LatLng_FromDegrees(extent[3], extent[2]))
        set_hex = [x.id() for x in coverer.GetCovering(region_rect)]

    else:
        coverer.set_fixed_level(res)
        region_rect = s2.S2LatLngRect(
        s2.S2LatLng_FromDegrees(-90, -180),
        s2.S2LatLng_FromDegrees(90, 180))
        set_hex = [x.id() for x in coverer.GetCovering(region_rect)]

    df = pd.DataFrame({"cell_id": set_hex})
    
    return df


def create_s2_geom_cells(extent, resolutions):


    # Create s2 rectangle to fill with s2 cells
    region_rect = s2.S2LatLngRect(
        s2.S2LatLng.FromDegrees(extent.bounds[1], extent.bounds[0]),
        s2.S2LatLng.FromDegrees(extent.bounds[3], extent.bounds[2]))
    coverer = s2.S2RegionCoverer()

    # Projection for cell area calculation
    transformer = Transformer.from_crs("epsg:4326", 'proj=isea')

    # Iterate through given resolutions, create and populate geopandas for each
    for r in resolutions:
        coverer.min_level = r
        coverer.max_level = r
        covering = coverer.GetCovering(region_rect)
        geoms = gpd.GeoDataFrame()
        geoms['cell_id'] = None
        geoms['area'] = None
        geoms['geometry'] = None
        for cellid in covering:
            new_cell = s2.S2Cell(cellid)
            vertices = []
            for i in range(0, 4):
                vertex = new_cell.GetS2LatLngVertex(i)
                vertices.append((vertex.lng().degrees,
                                 vertex.lat().degrees))
            geom = Polygon(vertices)
            geoms.loc[len(geoms)] = [cellid.get, transform(transformer.transform, geom).area, geom]

        geoms.to_file("s2_level{}.geojson".format(r), driver='GeoJSON')


def raster_to_s2(raster_path, value_name, cell_min_res, cell_max_res, extent=None, pix_size_factor=3):
    """Load raster values into s2 dggs cells

    Parameters:
    raster (string): path to raster file for uploading
    value_name (string): name of a value to be uploaded
    cell_min_res (integer): min h3 resolution to look for based on raster cell size
    cell_max_res (integer): max h3 resolution to look for based on raster cell size
    extent (list): Extent as array of 2 lon lat pairs to get raster values for
    pix_size_factor (pinteger): how times smaller h3 hex size should be comparing with raster cell size

    Returns:
    Pandas dataframe
   """
    # Open raster
    rs = rasterio.open(raster_path)

    # Get extent to fill with s2 squares
    if extent:
        region_rect = s2.S2LatLngRect(
            s2.S2LatLng_FromDegrees(extent[1], extent[0]),
            s2.S2LatLng_FromDegrees(extent[3], extent[2]))
    else:
        region_rect = s2.S2LatLngRect(
            s2.S2LatLng_FromDegrees(rs.bounds.bottom, rs.bounds.left),
            s2.S2LatLng_FromDegrees(rs.bounds.top, rs.bounds.right))

    # Get resolution dict
    resolutions = {}
    coverer = s2.S2RegionCoverer()
    # transformer = Transformer.from_crs("epsg:4326", 'proj=isea')
    for i in range(cell_min_res, cell_max_res, 1):
        # get s2 cell at level i
        coverer.set_fixed_level(i)
        cell = s2.S2Cell(coverer.GetCovering(region_rect)[0])

        # get s2 edge size at resolution i
        p1 = cell.GetS2LatLngVertex(0)
        p2 = cell.GetS2LatLngVertex(1)
        # edge = Point(transformer.transform(p2.lat().degrees(), p2.lng().degrees())).distance(Point(transformer.transform(p1.lat().degrees(), p2.lng().degrees())))
        edge = __haversine(p2.lat().degrees(), p2.lng().degrees(), p1.lat().degrees(), p2.lng().degrees())

        resolutions[i] = edge

    # Get two neighbour pixels in raster
    x1 = rs.transform[2]
    y1 = rs.transform[5]
    x2 = rs.transform[2] + rs.transform[0]
    y2 = rs.transform[5] - rs.transform[4]

    # Get pixel size from projected src
    size = __haversine(x1, y1, x1, y2)

    print(f"Raster pixel size {size}")

    # Get raster band as np array
    raster_band_array = rs.read(1)

    # Get h3 resolution for raster pixel size
    for key, value in resolutions.items():
        print(value)
        if value < size / pix_size_factor:
            resolution = key
            break
    print(resolution)

    coverer.set_fixed_level(resolution)

    # Create dataframe with cell_ids from cover with given resolution
    print(f"Start filling raster extent with s2 indexes at resolution {resolution}")
    df = pd.DataFrame({'cell_id': [x.id() for x in coverer.GetCovering(region_rect)]})

    # Get raster values for each hex_id
    print(f"Start getting raster values for s2 cells at resolution {resolution}")
    df[value_name] = df['cell_id'].apply(lambda x: raster_band_array[
        rs.index(s2.S2CellId(x).ToLatLng().lng().degrees(), s2.S2CellId(x).ToLatLng().lat().degrees())])

    # Drop nodata
    df = df[df[value_name] != rs.nodata]

    return df


def vector_to_s2(vector_path, value_name, resolution, extent=None, layer=None):
    """Load vector values into s2 dggs cells

    Parameters:
    vector_path (string): path to vector file for uploading
    value_name (string): name of a vector attribute to be uploaded
    resolution (integer): s2 resolution to load vector values into
    extent (list): Extent as array of 2 lat lon pairs to get vector values for

    Returns:
    Pandas dataframe
   """
    # Open vector to geodataframe
    gdf = gpd.read_file(vector_path, layer)

    # Get extent to fill with s2 squares
    if extent:
        region_rect = s2.S2LatLngRect(
            s2.S2LatLng_FromDegrees(extent[1], extent[0]),
            s2.S2LatLng_FromDegrees(extent[3], extent[2]))
    else:
        region_rect = s2.S2LatLngRect(
            s2.S2LatLng_FromDegrees(gdf['geometry'].total_bounds[1], gdf['geometry'].total_bounds[0]),
            s2.S2LatLng_FromDegrees(gdf['geometry'].total_bounds[2], gdf['geometry'].total_bounds[3]))

    coverer = s2.S2RegionCoverer()
    coverer.set_fixed_level(resolution)

    # Create dataframe with cell_ids from cover with given resolution
    print(f"Start filling raster extent with s2 indexes at resolution {resolution}")
    s2_gdf = gpd.GeoDataFrame({'cell_id': [x.id() for x in coverer.GetCovering(region_rect)]})

    # Get hex centroids for points
    s2_gdf['geometry'] = s2_gdf['cell_id'].apply(
        lambda x: Point(s2.S2CellId(x).ToLatLng().lng().degrees(), s2.S2CellId(x).ToLatLng().lat().degrees()))
    s2_gdf = s2_gdf.set_crs('epsg:4326')

    # Spatial join hex centroids with gdf
    vector_s2 = gpd.sjoin(s2_gdf, gdf)

    # Drop unnecessary fields
    vector_s2 = vector_s2[['cell_id', value_name]]

    return vector_s2


def cell_s2_downsampling(df, cell_id_col, metric_col, coarse_resolution, metric_type):
    """Aggregates a given attribute in s2 cell to a given coarser resolution level

    Parameters:
    df (pandas dataframe): dataframe with s2 ids and attributes for aggregation
    cell_id_col (string): name of s2 id column
    metric_col (string): name of a column for aggreagation
    coarse_resolution (integer): Coarser s2 resoluiton for aggregation
    metric_type (string): attribute type (numerical, categorical)
    Returns:
    Pandas dataframe
   """
    df_coarse = df.copy()
    coarse_id_col = 'cell_id_{}'.format(coarse_resolution)
    df_coarse[coarse_id_col] = df_coarse[cell_id_col].apply(lambda x: s2.S2CellId(x).parent(coarse_resolution).id())

    if metric_type == 'numeric':
        dfc = df_coarse.groupby(coarse_id_col)[[metric_col]].mean().reset_index()
    elif metric_type == 'categorical':
        dfc = df_coarse.groupby([coarse_id_col, metric_col]).agg(count=(metric_col, 'count')).reset_index().sort_values(
            by=[coarse_id_col, metric_col, 'count']).groupby(coarse_id_col, as_index=False, sort=False).first()
        dfc.drop('count', axis=1, inplace=True)

    dfc.columns = [cell_id_col, metric_col]

    return dfc