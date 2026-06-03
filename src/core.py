import re
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
from osgeo import gdal
import rasterio
import rasterio.warp
from rasterio import features
import numpy as np
import shutil

#####NOTEBOOK 1 DEM DOWNLOAD ###############

def nhd_download(select_data:str,vector_dir:pl.Path) -> str:
    '''
    creates link to nhd zip file on USGS The National Map. Checks link before attempting 
    to download. Prints updates and returns a list of the files downloaded.
    '''
    ##update in future to get S3 bucket information and only include available files
    s3_url = 'https://prd-tnm.s3.amazonaws.com'
    product_name = 'NHDPLUS_H_'+select_data[:4]+'_HU4_GDB.zip'
    link = s3_url+'/StagedProducts/Hydrography/NHDPlusHR/Beta/GDB/'+product_name
    gdb = product_name[:-4]+'.gdb'
    gdb_dir = vector_dir/gdb
    if vector_dir.exists() is False:
        os.mkdir(vector_dir)
    else:
        if gdb_dir.exists() is True:
            print('gdb is already downloaded')
            return gdb_dir
    print('opening: '+link)
    alive = url_is_alive(link)
    if alive is True:
        os.mkdir(gdb_dir)
        r = requests.get(link, allow_redirects=True)
        memfile = io.BytesIO(r.content)
        with ZipFile(memfile,'r') as openzip:
            print('saving to: '+str(gdb_dir))
            openzip.extractall(path=vector_dir)
        return gdb_dir
    else:
        print('download attempt failed')
        return 'attempt failed'

def url_is_alive(url:str)-> bool:
    '''
    Checks that a given URL is reachable.
    :param url: A URL
    :rtype: bool
    '''
    

    try:
        response = requests.get(url)
        if response.status_code == 200:
            return True
        else:
            return False
    except urllib.request.HTTPError:
        return False

def lower_pd_cols(df):
    columns = df.columns.to_list()
    lower_dict = {}
    for col in columns:
        lower_dict[col] = col.lower()
    df = df.rename(columns = lower_dict)
    return df

def tnm_coverage_check(gdf,aoi_coverage_buffer,projection,minimum_coverage = 0.95):
    'minimum % coverage at 95% hard coded'
    if gdf.empty:
        return False
    else:
        aoi_coverage = aoi_coverage_buffer.to_crs(epsg = projection)
        gdf_coverage = gdf.dissolve().to_crs(epsg = projection)
        inter_area = gdf_coverage.intersection(aoi_coverage)
        coverage_perc = inter_area.area.iloc[0]/aoi_coverage.area.iloc[0]
        if coverage_perc >= minimum_coverage:
            coverage = True
        else:
            coverage = False
        return coverage, coverage_perc
        

def tnm_coverage(aoi, desired_resolution='1m',buffer = 2000, minimum_coverage = 0.95):
    projection = '5070' #used for equal area across the US
    #create 2000' buffer
    aoi_coverage = aoi.dissolve().to_crs(epsg = projection)
    aoi_coverage_buffer = aoi_coverage.buffer(buffer).to_crs(epsg = '4269') #nad83
    if desired_resolution == 'opr':
        tiles=20
    else:
        tiles=4
    bounds = bound_box_str(aoi_coverage_buffer,tiles) #get bounding box  representation of geometry 
    DEM_products = {'1m':'Digital Elevation Model (DEM) 1 meter',
                    '3m': 'National Elevation Dataset (NED) 1/9 arc-second','10m':'National Elevation Dataset (NED) 1/3 arc-second Current', '30m':'National Elevation Dataset (NED) 1 arc-second',
                   'opr': 'Original Product Resolution (OPR) Digital Elevation Model (DEM)'}
    covered = False
    assert desired_resolution in list(DEM_products.keys()), "Desired DEM resolution can only be 1m, 3m, 10m, 30m, or opr" 
    for res in list(DEM_products.keys())[list(DEM_products.keys()).index(desired_resolution):]:
        tnm_responses = []
        dataset = DEM_products[res]
        for huc12_bound in bounds:
            tnm_out = tnm_query(dataset,huc12_bound)
            attmpt = 1
            while 'items' not in tnm_out.keys():
                print('missed api, retrying'+str(attmpt))
                attmpt += 1
                clear_output(wait=True)
                tnm_out = tnm_query(dataset,huc12_bound)
            if len(tnm_out['items']) == 100:
                print('need larger increment used. reset increment')
            tnm_responses += tnm_out['items']
        if tnm_responses:
            tnm_gdf_out = create_tnm_gdf(tnm_responses)
            tnm_gdf_out['dl_link'] = tnm_gdf_out['urls'].apply(lambda x: list(x.values())[0])
            tnm_gdf_clean = tnm_gdf_out.drop_duplicates(subset='dl_link')
            #trial by date
            collection_dates = list(tnm_gdf_clean["publicationDate"].unique())
            collection_dates.sort(reverse=True,key=lambda x: datetime.datetime.strptime(x, '%Y-%m-%d')) 
            #add dates until coverage check is met
            covered, coverage_perc = tnm_coverage_check(tnm_gdf_clean,aoi_coverage_buffer,projection,minimum_coverage)
            if covered is True:
                print('Best available DEM is '+res)
                print('coverage is %s percent' %(round(coverage_perc,3)*100))
                return tnm_gdf_clean, aoi_coverage_buffer
        else:
            covered, coverage_perc = [False,0.0]
        if covered is True:
            print('Best available DEM is '+res)
            print('coverage is %s percent' %(round(coverage_perc,3)*100))
            return tnm_gdf_clean, aoi_coverage_buffer
            
        else:
            print(res +' not available. Checking next best resolution')
            if res == 'opr':
                print('No high quality LiDAR coverage for this project area')
                return [None,None]

def bound_box_str(gdf:gpd.GeoDataFrame,increment:int):
    features = []
    ext_df = gdf.bounds.to_dict()
    idx = gdf.index.to_list()
    iterate = range(increment)
    for ix in idx:
        startY_big = ext_df['miny'][ix]
        startX_big = ext_df['minx'][ix]
        diffY = (ext_df['maxy'][ix] - startY_big)/increment
        diffX =  (ext_df['maxx'][ix] - startX_big)/increment
        for x,y in [(x,y) for x in iterate for y in iterate]:
            startY = startY_big + diffY*(y)
            startX = startX_big + diffX*(x)
            endY = startY_big + diffY*(y+1)
            endX = startX_big + diffX*(x+1)
            b_box = [str(round(startX,8)),str(round(endY,8)),str(round(endX,8)),str(round(startY,8))]
            features.append(', '.join(b_box))
    return features

def tnm_query(dataset:str,bounds:str) -> list:
    '''
    Uses TNM API service to get product download links 
    returns all relevant data from rest service and must
    be filtered by user to get desired content such as download URLs    
    '''
    import json
    base_url = 'https://tnmaccess.nationalmap.gov/api/v1/products?'
    
    attempts = 0
    while attempts < 50000:
        try:
            parameters = {'datasets':dataset,'bbox':bounds,'outputFormat':'JSON'}
            r = requests.get(base_url,parameters)
            response = json.loads(r.content)
            return response
        except:
            attempts+=1
            print(f'Good url but usgs service error. attempt {attempts}')
            print(r.url)
            bounds = bounds[:bounds.find(',')]+'0'+bounds[bounds.find(','):]
            clear_output(wait=True)

def create_tnm_gdf(responses:list):
    wkid = 4269 #NAD83
    epsg = 'EPSG:'+str(wkid)
    df = pd.DataFrame(responses)
    if df.empty:
        print('1m not available')
        return ['1m not available']
    else:
        b = df.apply(lambda row: geometry.box(row.boundingBox['minX'], row.boundingBox['minY'], row.boundingBox['maxX'], row.boundingBox['maxY']), axis=1)
        gdf = gpd.GeoDataFrame(df, geometry=b,crs = epsg)
    return gdf

def dem_download(dem_links:list,output_dir:str,dem_folder:list) -> list:
    '''
    Downloads links to DEM data. Does not check link before attempting 
    to download. Prints updates and returns a list of the dems downloaded.
    
    Need to add overwrite option
    '''
    dems = []
    total = len(dem_links)
    for link, fol in zip(dem_links,dem_folder):
        number = dem_links.index(link) + 1
        print('Obtaining DEM '+str(number)+ ' of '+str(total))
        output_dir_pl = pl.Path(output_dir)/fol
        if output_dir_pl.exists() is False:
            os.makedirs(output_dir_pl)

        raster_name = link.split('/')[-1]
        if (output_dir_pl/raster_name).exists() is False:
            print('opening: '+link)
            r = requests.get(link, allow_redirects=True)
            print('saving to: '+str(output_dir_pl/raster_name))
            with open (output_dir_pl/raster_name, 'wb') as tif:
                tif.write(r.content)
            dems.append(raster_name)
           
        else:
            print('DEM already downloaded')
            dems.append(raster_name)
        clear_output(wait=True)
    return dems

def process_dems(raster_paths, select_data, sr, overwrite=False):
    """
    merges all tif rasters to VRT and then reprojects to the given spatial reference.
    """
    raster_dir = pl.Path(raster_paths[0]).parent.parent
    out_vrt_big = pl.Path('{0}/dem_{1}.vrt'.format(raster_dir,select_data))
    for raster in raster_paths:
        assert os.path.exists(raster), "raster doesn't exist: " + str(raster) 
    
    if overwrite == False and check_files_exist([out_vrt]): 
        print("processed DEM already exists; not overwriting.")
        return 
    
    pl.Path(raster_dir.parent).mkdir(parents=True, exist_ok=True)
    first_raster = raster_paths[0]
    if not sr:
        sr = get_raster_epsg(first_raster)
    temp_name = str(pl.Path(raster_dir).parent/str(first_raster).split("/")[-1])
    ext = temp_name[temp_name.find('.'):]
    for raster in raster_paths:
        name = str(raster)
        reproject_raster(
            raster, 
            out_vrt = name.replace(ext, ".vrt"),
            sr=sr, 
            overwrite=overwrite
        )

    reprojected_rasters = [str(raster).replace(ext, ".vrt") for raster in raster_paths]     
    merged_rasters = gdal.BuildVRT(str(out_vrt_big), reprojected_rasters)
    return out_vrt_big

def get_raster_epsg(file_path):
    with rasterio.open(file_path) as src:
        sr = str(src.crs.to_epsg())
    return sr

def reproject_raster(raster_path:list, out_vrt, sr:str, overwrite=False):
    '''
    Reprojects rasters to provided spatial reference and saves output as vrt file. 
    Returns a list of the paths to the reporjected rasters
    
    :param raster_paths: full paths to files to re-project
    :param str sr: EPSG spatial reference to reproject to
    :param str out_dir: output path to save rasters to
    '''
    assert os.path.exists(raster_path), "raster doesn't exist missing: " + str(raster_path)
    if overwrite is False and check_files_exist([out_vrt]): 
        print("files already exists; not recalculating.")
        return
    if len(sr) > 4: 
        to_crs = pyproj.CRS.from_authority("esri", int(sr)).to_wkt()
        pass
    else:
        to_crs = 'EPSG:'+ sr
    kwargs = gdal.WarpOptions(dstSRS =to_crs, resampleAlg= 'bilinear', errorThreshold=0.01)
    gdal.Warp(
        str(out_vrt), 
        str(raster_path),
        options= kwargs)

def check_dem_coverage(dem_path, huc_path,minimum_coverage = 0.95):
    """check if the DEM has full coverage of the HUC by 
    confirming the HUC doesn't intersect the DEM bounding box. 
    returns True if DEM covers full HUC, False otherwise"""
    projection = '5070' #equal area used for coverage check
    
    #get dem coverage approximation
    out_chk = dem_path.parent/'{0}_coverage.tif'.format(dem_path.name[:-4])
    with rasterio.open(dem_path) as src:
        kwds = src.profile
        warped = rasterio.warp.aligned_target(src.profile['transform'], src.profile['width'], src.profile['height'], src.profile['transform'][0]*src.profile['blockxsize'])
        kwds['driver'] = 'GTiff'
        kwds['dtype'] = rasterio.float32
        kwds['width'] = warped[1]
        kwds['height'] = warped[2]
        kwds['transform'] = warped[0]
        kwds['nodata'] = np.nan
        null_val = src.profile['nodata']
        coverage = np.array([[2]*kwds['width']]*kwds['height'],dtype=np.float32)
        for ji, window in src.block_windows(1):
            r = src.read(1, window=window)
            #slice by 10ths
            small = r[::10,::10]
            #test to see if all sliced values are null
            if np.nanmax(small) == null_val:
                coverage[ji[0],ji[1]] = np.nan
            else:
                coverage[ji[0],ji[1]] = 1
        with rasterio.open(out_chk, 'w', **kwds) as wrt:
            wrt.write(coverage,1)

    gpd_polygonized_raster = polygonize_raster(out_chk,projection)
    dem_bounds = gpd_polygonized_raster.loc[gpd_polygonized_raster['raster_val']== 1].dissolve()
    dem_bounds_normalized = dem_bounds.buffer(0)
    out_chk_vect = out_chk.parent/'{}.geojson'.format(out_chk.name[:-4])
    dem_bounds_normalized.to_file(out_chk_vect, driver="GeoJSON")
    os.unlink(out_chk)
    
    # read huc boundary 
    huc_geo = gpd.read_file(huc_path)
    huc_geo_sr = huc_geo.to_crs(dem_bounds_normalized.crs)
    huc_geo_sr.reset_index(inplace=True)

    #plot
    fix1, ax1 = plt.subplots(figsize=(30,10))

    huc_geo_sr.plot(ax = ax1, color = 'orange', edgecolor = 'white')
    dem_bounds_normalized.plot(ax=ax1, color='red',edgecolor = 'black',alpha=0.3)
    plt.show()

    #check coverage
    inter_area = dem_bounds_normalized.intersection(huc_geo_sr)
    coverage_perc = inter_area.area.iloc[0]/huc_geo_sr.area.iloc[0]
    if coverage_perc >= minimum_coverage:
        coverage = True
    else:
        coverage = False
    print('coverage is %s percent' %(round(coverage_perc,3)*100))
    return coverage

def polygonize_raster(raster_path,projection):
    mask = None
    with rasterio.Env():
        with rasterio.open(raster_path) as src:
            origin_sr = src.crs
            image = src.read(1) # first band
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                features.shapes(image, mask=mask, transform=src.transform)))
           
    geoms = list(results)
    output_gdf = gpd.GeoDataFrame.from_features(geoms,crs=origin_sr)
    #return output_gdf.to_crs(epsg=projection)
    return output_gdf