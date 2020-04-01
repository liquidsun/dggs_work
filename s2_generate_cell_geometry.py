import pyproj
from s2sphere import *
from shapely.geometry import Polygon
from shapely.ops import transform
import geopandas as gpd
from fiona.crs import from_epsg
from functools import partial

def create_s2_geom_cells(extent, resolutions):
    # Create s2 rectangle to fill with s2 cells
    region_rect = LatLngRect(
        LatLng.from_degrees(extent.bounds[1], extent.bounds[0]),
        LatLng.from_degrees(extent.bounds[3], extent.bounds[2]))
    coverer = RegionCoverer()

    # Projection for cell area calculation
    proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'),
                   pyproj.Proj('proj=isea'))

    # Iterate through given resolutions, create and populate geopandas for each
    for r in resolutions:
        coverer.min_level = r
        coverer.max_level = r
        covering = coverer.get_covering(region_rect)
        geoms = gpd.GeoDataFrame()
        geoms['s2id'] = None
        geoms['area'] = None
        geoms['geometry'] = None
        for cellid in covering:
            new_cell = Cell(cellid)
            vertices = []
            for i in range(0, 4):
                vertex = new_cell.get_vertex(i)
                latlng = LatLng.from_point(vertex)
                vertices.append((latlng.lng().degrees,
                                 latlng.lat().degrees))

            geom = Polygon(vertices)
            geoms.loc[len(geoms)] = [cellid.to_token(), transform(proj, geom).area, geom]

        geoms.to_file("s2_level{}.geojson".format(r), driver='GeoJSON')



