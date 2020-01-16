#!/usr/bin/env python

# T. Carman Winter 2019-2020

# Interactive viewer for mapping collections of dvm-dos-tem input datasets.


import gdal  # does not work on vagrant, or osx python 2.7x, works on atlas, and modex I think and OSX pyhton 3.7
import osr
import os

import geojson
import geopandas
import json
import glob

from bokeh.io import output_file, show, curdoc
from bokeh.models.sources import GeoJSONDataSource
from bokeh.models import HoverTool, TapTool, Patches
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure

from bokeh.tile_providers import get_provider, Vendors


def get_corner_coords(file_path):
  '''
  Given a file path to a vegetation.nc dvmdostem input file, this function
  figured out the corner coordinates in WGS84. Assumes that the vegetation.nc
  file is projected and has all the spatial reference info.
  Returns coordinate pairs for lr, ll, ul, ur.
  '''
  ds = gdal.Open("NETCDF:{}:veg_class".format(file_path))
  ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()

  lrx = ulx + ds.RasterXSize * xres
  lry = uly + ds.RasterYSize * yres

  llx = ulx
  lly = uly + ds.RasterYSize * yres

  urx = ulx + ds.RasterXSize * xres
  ury = uly

  targetSR = osr.SpatialReference()
  targetSR.ImportFromEPSG(4326) # WGS84
  sourceSR = osr.SpatialReference(wkt=ds.GetProjectionRef())
  coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)

  return [coordTrans.TransformPoint(X, Y) for X,Y in zip((lrx,llx,ulx,urx),(lry,lly,uly,ury))]


def get_io_folder_listing(starting_path, pattern):
  '''Get list of folders using the starting path and matching the pattern.'''
  folder_list = glob.glob(os.path.join(starting_path, pattern))
  folder_list = [i for i in folder_list if ".tar.gz" not in i]
  folder_list = [i for i in folder_list if os.path.isdir(i)]

  return folder_list


def site_table(folder_list):
  '''
  Print table with corner coords for each site. Assumes that each folder in 
  `folder_list` contains a vegetation.nc file which is geo-referenced.

  Table looks like this:
    site,lr_lon,lr_lat,ll_lon,ll_lat,ul_lon,ul_lat,ur_lon,ur_lat
    cru-ts40_ar5_rcp85_ncar-ccsm4_CALM_Westdock8_Met_10x10,-148.3557,70.3419,-148.6162,70.3495,-148.5949,70.4406,-148.3335,70.4330
    cru-ts40_ar5_rcp85_ncar-ccsm4_SNOTEL_SUSITNAVALLEYHIGH_10x10,-149.8489,62.1247,-150.0409,62.1302,-150.0294,62.2194,-149.8369,62.2140
  '''
  print("site,lr_lon,lr_lat,ll_lon,ll_lat,ul_lon,ul_lat,ur_lon,ur_lat")
  for f in folder_list:
    f = os.path.join(f, 'vegetation.nc')
    
    lr,ll,ul,ur = get_corner_coords(f)
    print("{},{:0.4f},{:0.4f},{:0.4f},{:0.4f},{:0.4f},{:0.4f},{:0.4f},{:0.4f}".format(
        os.path.basename(os.path.dirname(f)),
        lr[0], lr[1],
        ll[0], ll[1],
        ul[0], ul[1],
        ur[0], ur[1],
      )
    )

def build_feature_collection(folder_list):
  '''
  Builds a GeoJSON feature collection with a polygon for site in the `folder_list`.
  '''

  geojson_obj = dict(type="FeatureCollection",features=[])

  for folder in folder_list:
    f = os.path.join(folder, 'vegetation.nc')
    
    lr,ll,ul,ur = get_corner_coords(f)

    geojson_obj['features'].append(
      dict(
        type="Feature", 
        geometry=dict(
          type="Polygon", 
          coordinates=[
            [list(lr),list(ll),list(ul),list(ur),list(lr)]
          ]
        ), 
        properties=dict(
          name=f
        )
      )
    )

  return geojson_obj


def build_geojson(feature_collection, method='A', proj='epsg:3857'):
  if method == 'A':
    # Write sample data to file
    ff = "/tmp/junk.geojson"
    with open(ff, 'w') as f:
      f.write(geojson.dumps(feature_collection))

    # Read data from file
    gj = geopandas.read_file(ff)
    os.remove("/tmp/junk.geojson")

  if method == 'B':
    # This messes up something with the CRS...doesn't get created correctly.
    # Says we have to set it first...
    gj = geopandas.GeoDataFrame.from_features(feature_collection)

  # Re-project data
  gj = gj.to_crs(epsg=3857) # Convert to web mercator projection that tile provider uses
  return gj

#for i in range(len(data["features"])):
#  data['features'][i]["properties"]["Color"] = ["blue","red"][i%2]

def tapselect_handler(attr, old, new):
  print("attr:{}   old:{}    new:{}".format(attr, old, new))
  print("Somewhere to stop!")
  print(input_areas.data_source)
  from IPython import embed; embed()
  #import pdb; pdb.set_trace()

#### Prepare data

folders = get_io_folder_listing("/Users/tobeycarman/Documents/SEL/dvmdostem-input-catalog", pattern="*")

feature_collection = build_feature_collection(folders)

my_geojson = build_geojson(feature_collection) 


# Discover bounds
xmin, ymin, xmax,ymax=my_geojson.total_bounds

# Conver to Bokeh friendly object
geosource = GeoJSONDataSource(geojson=my_geojson.to_json())


##### Plot
output_file("geojson.html")

p = figure(x_range=(xmin, xmax), y_range=(ymin, ymax),
           x_axis_type="mercator", y_axis_type="mercator")

# Add tiles
tile_provider = get_provider(Vendors.CARTODBPOSITRON)
p.add_tile(tile_provider)


# Have to add the patches after the tiles
input_areas = p.patches('xs','ys', source=geosource, fill_color='Red', line_color='gray', line_width=0.25, fill_alpha=.3)

taptool = TapTool(renderers = [input_areas])
hovertool = HoverTool(renderers = [input_areas],
                      tooltips = [
                        ('Site','@name')
                      ])

selected_patches = Patches(fill_color="yellow", fill_alpha=.8, line_width=1, line_color="green")
input_areas.selection_glyph = selected_patches

geosource.selected.on_change('indices', tapselect_handler)


# CustomJS(code="""
# // comment....
# console.log("Somewhere to stop!");
# """)

# Create hover tool
p.add_tools(taptool, hovertool)

curdoc().add_root(p)



#show(p)








