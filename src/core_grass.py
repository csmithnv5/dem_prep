import re
import sys
import os
import pathlib as pl
import requests
import io
from zipfile import ZipFile
import pandas as pd
import geopandas as gpd
import pyproj
from shapely import geometry
from shapely.geometry import shape, Polygon
import datetime
import glob
from IPython.display import clear_output
import numpy as np
import shutil

import subprocess
result = subprocess.run(["grass", "--config", "python_path"],
                        check=True, text=True, capture_output=True)
sys.path.append(result.stdout.strip())

# import GRASS Python packages
from collections import defaultdict
import grass.script as gs
import grass.jupyter as gj

###NOTEBOOK 2 ######
def list_existing_grass(print_it=False):
    #List Existing Files: Vectors and Rasters
    layers = {'vector':[],'raster':[]}
    if print_it:
        print('Available vector maps:')
    for vect in gs.list_strings(type='vector'):
        if print_it:
            print (vect)
        layers['vector'].append(vect)

    if print_it:
        print('\nAvailable raster maps:')
    for rast in gs.list_strings(type='raster'):
        if print_it:
            print (rast)
        layers['raster'].append(rast)
    return layers

def dem_tiles_to_gis(gs,dem_dir,Project_Area,vert_unit, patch=False,delete_raw=False,force=False):
    """imports rasters to given grass session in gs"""
    GRASS_raster_files= [file for line in gs.list_strings(type='raster') for file in [line.split("@")[0]]]
    dem_name = Project_Area + '_dem'
    dem_boundary_v = dem_name+'_boundary'
    if dem_name in GRASS_raster_files:
        print('dem already created, skipping import and patching. Use Force = True to recreate dem')
        return dem_name
    else:
        tile_rasters = []
        dem_b_rasters = [] 
        clean_up_rasters = [] #hold raster list for later cleanup
        li_all_files = raster_files_dic(dem_dir)
        if not li_all_files:
            print('dem already imported and deleted. Must redownload')
            return None
        total_rasters = sorted({x for v in li_all_files.values() for x in v})
        for collection, rasters in li_all_files.items():
            for raster in rasters:
                name = ('%s' % pl.Path(raster).name)
                ext = name[name.find('.'):]
                raster_name = ('%s' % pl.Path(raster).name)[:-4]
                vrt_path = os.path.join(pl.Path(raster).parent, raster_name + ".vrt")
                raster_name_u = f'{raster_name}_{collection}'.replace('-','_')
                if raster_name not in GRASS_raster_files:
                    print("adding raster {}".format(raster_name_u))
                    gs.run_command('r.in.gdal', input=vrt_path, output= raster_name_u)
                    gs.run_command('g.region', raster=raster_name_u) 
                    dem_vert_unit_check(raster,raster_name_u,vert_unit)

                tile_rasters.append(raster_name_u)
                if delete_raw is True:
                    os.unlink(raster)

                clear_output(wait=True)
        if not tile_rasters:
            for dir_path, dir_names, file_names in os.walk(dem_dir):
                tile_rasters += [file.split(ext)[0] for file in file_names]
        if patch is True:
            print("Patching together {} rasters".format(len(tile_rasters)))

            if len(total_rasters) > 1:
                gs.run_command('g.region', raster=tile_rasters)
                gs.run_command('r.patch', input=tile_rasters, output=dem_name)

            else:
                gs.run_command('g.rename', raster=(tile_rasters[0],dem_name))
    #             gs.run_command('g.rename', raster=(dem_b_rasters[0],dem_name))

            ##fill any gaps in the dem
            #gs.run_command('r.fill.gaps',input=DEM,output='tmp', method)

            print('Created: '+ dem_name)
    #         print('Created: '+ dem_name+'_b')

            #cleanup and remove raster tiles
            for tile in tile_rasters:
                gs.run_command('g.remove', type ='raster', name=tile, flags ='f')
    #             gs.run_command('g.remove', type ='raster', name=tile+'_b', flags ='f')
            return dem_name
        else:
            return tile_rasters

def raster_files_dic(path: str):
    '''Create a list in the directory given by the path string. Looks for all tiff and img files. 
    adds them in order of date
    '''
    li_raster_files = defaultdict(list)
    folders = [name for name in os.listdir(path) if os.path.isdir(path/name)]
    
    #remove hidden folders on mac file systems
    for f in folders:
        if f.find('.')>=0:
            folders.remove(f)
            
    folders.sort(key=lambda f: int(re.sub('\D', '', f)), reverse = True)
    for f_dir in folders:
        li_files = [file for file in os.listdir(path/f_dir) if os.path.splitext(file)[-1] in ['.tif','.img']]
        for file in li_files:
            li_raster_files[f_dir].append(os.path.join(path,f_dir, file))
    return li_raster_files

def dem_vert_unit_check(dem_source,raster_name,vert_unit):
    vert_value = {'feet':float(0.3048),'meter':float(1)}
    meta = gs.parse_command('g.proj', georef= dem_source, flags='g')
    vert_value_m = float(meta['meters'])
    vert_conversion = (vert_value_m/vert_value[vert_unit])
    if vert_conversion != float(1):
        print('converting vertical units to '+vert_unit+'. Raw DEM * '+str(vert_conversion))
        gs.run_command('g.region', raster=raster_name)
        gs.mapcalc("new_rast = {0} * {1}".format(raster_name, vert_conversion))
        gs.run_command('g.rename', raster= "new_rast,{}".format(raster_name))
        gs.run_command('r.support', map= raster_name, units='ft')
    else:
        gs.run_command('r.support', map= raster_name, units='m')
        
    return None

def basin_dems(DEM_tiles,basins,basin_list,basin_scale,raster_dir,buffer = 0,fill_nulls = True, remove_tiles = False,export_tiff = True,force=True):
    basin_DEMs = []
    for basin in basin_list:
        dem_name = 'dem_'+str(basin)
        gtiff = dem_name+'.tif'
        
        if gtiff not in os.listdir(raster_dir) or force:
            if type(basin) == str:
                gs.run_command('v.extract', input =basins, where= "%s = '%s'" % (basin_scale,basin), output = 'tmp_bsn')
            else:
                gs.run_command('v.extract', input =basins, where= "%s = %s" % (basin_scale,basin), output = 'tmp_bsn')
            
            gs.run_command('v.buffer',input='tmp_bsn',output='tmp_bsn_buffer',distance = buffer)
               
            print("Patching together rasters for", basin)
            #set the region based on the buffer but align with DEM raster
            gs.run_command('g.region',vector = 'tmp_bsn_buffer',align=DEM_tiles[0])
            #mask outside cells
            gs.run_command('r.mask',vector = 'tmp_bsn_buffer')
            
            if len(DEM_tiles) > 1:
                gs.run_command('r.patch', input=DEM_tiles, output=dem_name)

            else:
                gs.run_command('g.rename', raster=(DEM_tiles[0],dem_name))
            
            if fill_nulls == True:
                gs.run_command('r.fillnulls', input = dem_name,output = 'dem_tmp', method = 'bilinear')
                gs.run_command('g.rename', raster=('dem_tmp',dem_name))

            #export to gdal
            if export_tiff == True:
                gs.run_command('r.out.gdal', input= dem_name, output= raster_dir/gtiff,createopt="BIGTIFF=YES,PROFILE=GeoTIFF,TFW=YES,COMPRESS=LZW",flags='f')
        
        print('created '+gtiff)
        basin_DEMs.append(dem_name)
    
    GRASS_raster_files= [file for line in gs.list_strings(type='raster') for file in [line.split("@")[0]]]
    
    #remove intermediate tiles option
    if remove_tiles:
        for tile in DEM_tiles:
            gs.run_command('g.remove', type ='raster', name=tile, flags ='f')
    if 'MASK' in GRASS_raster_files:
        gs.run_command('r.mask',flags='r')
    return basin_DEMs