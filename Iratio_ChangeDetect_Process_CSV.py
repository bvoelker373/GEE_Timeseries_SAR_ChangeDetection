# Initialize Earth Engine Python API.
import ee
#ee.Authenticate()
ee.Initialize()

# For data manipulation and analysis.
import math
import pandas as pd
import numpy as np
np.set_printoptions(precision=4, suppress=True)
import scipy.signal

# For plotting
import matplotlib.pyplot as plt
from matplotlib import gridspec

def MakeGrid(geometry, scale):
    """Takes a polygon and creates a grid of polygons inside the shape.
    Adapted from: https://developers.google.com/earth-engine/tutorials/community/drawing-tools
    
    Keyword arguments:
    geometry -- Earth Engine geometry object
    scale -- desired spacing of grid
    """
    # pixelLonLat returns an image with each pixel labeled with longitude and
    # latitude values.
    lonLat = ee.Image.pixelLonLat()

    # Select the longitude and latitude bands, multiply by a large number then
    # truncate them to integers.
    lonGrid = lonLat.select('longitude').multiply(10000000).toInt()

    latGrid = lonLat.select('latitude').multiply(10000000).toInt()

    # To produce the grid, multiply the latitude and longitude images and then use
    # reduce to vectors at the 10km resolution to group the grid into vectors.
    return lonGrid.multiply(latGrid).reduceToVectors(geometry = geometry, scale = scale, geometryType = 'polygon',)

def Add_ID_to_Features(dataset):
    """Gives a unique ID number to each feature in a feature collection, as a new property.
    Adapted from: https://gis.stackexchange.com/questions/374137/add-an-incremential-number-to-each-feature-in-a-featurecollection-in-gee
    
    Keyword argument:
    dataset -- Earth Engine feature collection
    """
    indexes = ee.List(dataset.aggregate_array('system:index'))
    feat_ids = ee.List.sequence(1, indexes.size())
    idByIndex = ee.Dictionary.fromLists(indexes, feat_ids)
    
    def function(feature):
        """Adds the ID number to the feature."""
        return feature.set('ID', idByIndex.get(feature.get('system:index')))
    # Map the function over the dataset.
    return dataset.map(function)

def Get_BeforeAfter_Imagery(out_dir, ID, event_time, region, sizeWindows):
    """Obtain cloud-masked Sentinel-2 imagery before/after a given date.
    
    Keyword arguments:
    ID -- unique integer identifier of polygon representing landslide detection area
    event_time -- mean date of detection window (input format)
    region -- bounding box of area of interest (input format)
    sizeWindows -- duration of the detection window in days
    """
    def MaskS2clouds(image):
        """Filter and mask clouds for Sentinel-2 optical data
        
        From: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2
        """
        qa = image.select('QA60')
        #Bits 10 and 11 are clouds and cirrus, respectively.
        cloudBitMask = 1 << 10
        cirrusBitMask = 1 << 11
        # Both flags should be set to zero, indicating clear conditions.
        mask = (qa.bitwiseAnd(cloudBitMask).eq(0)
                .And(qa.bitwiseAnd(cirrusBitMask).eq(0))
               )
        return image.updateMask(mask).divide(10000)

    str_pre = "{}_{}_pre".format(ID, str(event_time)[ 0 : 10 ])
    str_post = "{}_{}_post".format(ID, str(event_time)[ 0 : 10 ])
    
    event_time_ee = ee.Date(event_time)
    
    # Pre-window time period for pre-event S2 image collection.
    preWindow_T1 = event_time_ee.advance(-sizeWindows - 30, 'day')
    preWindow_T2 = event_time_ee.advance(-sizeWindows, 'day')

    # Post-window time period for post-event S2 image collection.
    postWindow_T1 = event_time_ee.advance(sizeWindows, 'day')
    postWindow_T2 = event_time_ee.advance(sizeWindows + 30, 'day')

    optical_pre = (ee.ImageCollection('COPERNICUS/S2')
                   .filterDate(preWindow_T1, preWindow_T2)
                   .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
                   .filterBounds(region)
                   .map(MaskS2clouds)
                   .median()
                  )

    optical_post = (ee.ImageCollection('COPERNICUS/S2')
                    .filterDate(postWindow_T1, postWindow_T2)
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
                    .filterBounds(region)
                    .map(MaskS2clouds)
                    .median()
                   )
    
    exportRegion = ee.Geometry.Polygon(region.filter(ee.Filter.eq('ID', ID)).first().getInfo()['geometry']['coordinates'])
    
    exportImgPre = ee.batch.Export.image.toDrive(image = optical_pre,
                                                 folder = out_dir,
                                                 description = str_pre,
                                                 region = exportRegion,
                                                 scale = 10,
                                                 maxPixels = 1e9,
                                                 fileFormat = 'GeoTIFF')
    exportImgPost = ee.batch.Export.image.toDrive(image = optical_post,
                                                  folder = out_dir,
                                                  description = str_post,
                                                  region = exportRegion,
                                                  scale = 10,
                                                  maxPixels = 1e9,
                                                  fileFormat = 'GeoTIFF')
    exportImgPre.start()
    exportImgPost.start()

def Get_BeforeAfter_NDVI(ID, event_time, region, sizeWindows):
    """Obtain cloud-masked median NDVI before/after a given date.
    
    Keyword arguments:
    ID -- unique integer identifier of polygon representing landslide detection area
    event_time -- mean date of detection window (input format)
    region -- bounding box of area of interest (input format)
    sizeWindows -- duration of the detection window in days
    """
    def AddNDVI(image):
        """Adds an NDVI band to a Sentinel-2 image. NDVI is calculated as
        the normalized difference between the near-infrared band and the red band,
        which correspond to the 8th and 4th band in the Sentinel-2 imagery.
        
        From: https://developers.google.com/earth-engine/tutorials/tutorial_api_06
        """
        ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
        return image.addBands(ndvi)

    # Function to filter and mask clouds for Sentinel-2 optical data
    def MaskS2clouds(image):
        """Filter and mask clouds for Sentinel-2 optical data
        
        From: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2
        """
        qa = image.select('QA60')
        #Bits 10 and 11 are clouds and cirrus, respectively.
        cloudBitMask = 1 << 10
        cirrusBitMask = 1 << 11
        # Both flags should be set to zero, indicating clear conditions.
        mask = (qa.bitwiseAnd(cloudBitMask).eq(0)
                .And(qa.bitwiseAnd(cirrusBitMask).eq(0))
               )
        return image.updateMask(mask).divide(10000)
    
    event_time_ee = ee.Date(event_time)
    
    # Define pre-event time period, for pre-event NDVI image.
    preWindow_T1 = event_time_ee.advance(-sizeWindows - 30, 'day')
    preWindow_T2 = event_time_ee.advance(-sizeWindows, 'day')

    # Define post-event time period, for post-event NDVI image
    postWindow_T1 = event_time_ee.advance(sizeWindows, 'day')
    postWindow_T2 = event_time_ee.advance(sizeWindows + 30, 'day')

    # Get Sentinel 2 surface reflectance imagery before and after the window.
    s2_sr_before = ee.ImageCollection('COPERNICUS/S2').filterDate(preWindow_T1, preWindow_T2).filterBounds(region)
    s2_sr_after = ee.ImageCollection('COPERNICUS/S2').filterDate(postWindow_T1, postWindow_T2).filterBounds(region)

    # Apply the cloud masking function to the before/after image collections.
    s2_sr_before = s2_sr_before.map(MaskS2clouds)
    s2_sr_after = s2_sr_after.map(MaskS2clouds)

    # Apply the NDVI function to the before/after image collections.
    s2_ndvi_before = s2_sr_before.map(AddNDVI)
    s2_ndvi_after = s2_sr_after.map(AddNDVI)

    # Find median of the images in the pre-event image collection to get a pre-event image.
    pre_event_NDVI_img = s2_ndvi_before.select('NDVI').median()
    # Find median of the images in the post-event image collection, to get a post-event image.
    post_event_NDVI_img = s2_ndvi_after.select('NDVI').median()
    
    # Get the average NDVI over the area
    pre_NDVI = pre_event_NDVI_img.reduceRegion(
                                   geometry = region.filter(ee.Filter.eq('ID', ID)),
                                   reducer = ee.Reducer.mean(),
                                   scale = 10,
                                   bestEffort = True,
                                   maxPixels = 1e9
                                   )
    
    post_NDVI = post_event_NDVI_img.reduceRegion(
                                     geometry = region.filter(ee.Filter.eq('ID', ID)),
                                     reducer = ee.Reducer.mean(),
                                     scale = 10,
                                     bestEffort = True,
                                     maxPixels = 1e9
                                     )
    
    return pre_NDVI, post_NDVI


# Output Drive folder for CSVs and GeoTiffs.
out_dir = 'EE_SAR_MovingWindow_ChangeDetection_HolyFire'

# Define duration of window.
sizeWindows = 7 # days

# Area of interest to make a grid over.
ee_import = ee.FeatureCollection("users/bvoelker/Wildfires/Holy_Fire_2018_08_06_Boundary")
#ee_import = ee.FeatureCollection("users/bvoelker/hiroshima_landslide_subset2")
aoi_bounding_box = ee_import.geometry().bounds();

# The polygons to check for changes in will be a regularly spaced grid.
grid_size = 500 # 500m scale
grid = MakeGrid(aoi_bounding_box, grid_size) 
grid = Add_ID_to_Features(grid)

# Number of cells in grid.
numPolys = grid.size().getInfo()

# Load in the CSVs:
# Contains many irrelevant columns. The needed columns are the shapefile
#   joinkeys and the intensity ratio data of each polygon.
data_sum_raw = pd.read_csv('HF_CA_99th_Ptile_sum.csv')
data_count_raw = pd.read_csv('HF_CA_pixel_count.csv')

# Read in window dates
windows = pd.read_csv('windows_7day.csv')
windows['Time'] = pd.to_datetime(windows['Time'])

# Create an array of each polygon's unique identifier.
ids = data_sum_raw['ID'][0:numPolys]

# Format the mean intensity ratio:
# Rows: Zones (Landslide polygons).
# Columns: Mean Intensity ratio at each moving window date.
data_sum = np.transpose(data_sum_raw['sum']
                        .values
                        .reshape(-1, numPolys)
                       )
data_count = np.transpose(data_count_raw['count']
                          .values
                          .reshape(-1, numPolys)
                         )

detection_date_list = []
detection_ID_list = []
detection_index_list = []


# Dividing by a count of 0 results in NaN and a warning. Supress the warning.
np.seterr(all="ignore")
data_normalized = data_sum / data_count
# Change NaN values to 0.
data_normalized = np.nan_to_num(data_normalized)

for row_number, row in enumerate(data_normalized):
    """Adapted from: https://pythonawesome.com/overview-of-the-peaks-dectection-algorithms-available-in-python/"""
    peak_index = scipy.signal.argrelextrema(row,
                                            comparator = np.greater,
                                            order = 3)
    # Get first date where:
    # 1) Peak is detected
    # 2) Above pixel percentage threshold of 10%
    # 3) 7 months have passed since last detection. (?)
    detect = False
    detect_index = None
    for index in peak_index[0]:
        if detect is False:
            if row[index] > 0.2:
                detect = True
                detect_index = index
                detection_date_list.append(windows['Time'][detect_index])
                detection_ID_list.append(int(ids[row_number]))

# Get list of indexes too, which start at 0 (offset of 1 from the IDs)
detection_index_list = [x - 1 for x in detection_ID_list]


"""Change the figsize tuple to suit the plot as needed."""
figsize = (30,20)

numDetections = len(detection_ID_list)

plt_rows = math.ceil(math.sqrt(numDetections))
plt_cols = math.floor(math.sqrt(numDetections))

fig, axes = plt.subplots(plt_rows, plt_cols, sharex=True, sharey=True, figsize = figsize)

for i, ax in enumerate(axes.flatten()):
    if i > numDetections - 1:
        pass
    else:
        currentData = data_normalized.T[:,detection_index_list[i]]
        newLabel = "Cell #" + str(detection_ID_list[i])
        ax.plot(windows['Time'], currentData, color='black', label = newLabel)
        ax.axhline(0.2, alpha = 0.2, color = 'k', linestyle = '--')
        
        """Adapted from: https://pythonawesome.com/overview-of-the-peaks-dectection-algorithms-available-in-python/"""
        peak_index = scipy.signal.argrelextrema(currentData,
                                                comparator = np.greater,
                                                order = 3)
        
        ax.plot(windows['Time'][peak_index[0]], currentData[peak_index[0]], 'ro')
        ax.legend()
# Save time series figure.
plt.savefig('normalized_99Ptl_timeseries_HF.png')

plt.show()


# Get an image of every cell with a detected change before and after the detection date.
"""Will begin export of (potentially on the order of 10-100) images to Google Drive
Comment out the line below if export is not desired."""
[Get_BeforeAfter_Imagery(ID, date, grid, sizeWindows) for ID, date in zip(detection_ID_list, detection_date_list)]

# Get NDVI of every cell with a detected change before and after the detection date.
pre_NDVI = []
post_NDVI = []
for detection in range(len(detection_ID_list)):
    detID = detection_ID_list[detection]
    detDate = detection_date_list[detection]
    
    new_pre_NDVI, new_post_NDVI = Get_BeforeAfter_NDVI(detID, detDate, grid, sizeWindows)
    pre_NDVI.append(new_pre_NDVI.getInfo()['NDVI'])
    post_NDVI.append(new_post_NDVI.getInfo()['NDVI'])

NDVI_data = np.column_stack((pre_NDVI, post_NDVI))
np.savetxt("NDVI_changes.csv", NDVI_data, delimiter=",")
print(NDVI_data)
