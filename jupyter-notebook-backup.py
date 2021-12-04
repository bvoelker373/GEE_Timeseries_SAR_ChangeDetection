# To access Earth Engine Python API.
import ee
ee.Authenticate()
ee.Initialize()

# For data manipulation and analysis.
import math
import pandas as pd
import numpy as np
np.set_printoptions(precision=4, suppress=True)
from datetime import datetime
import scipy.signal

# For plotting
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.dates as mdates
import datetime as dt


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

def WindowDate_to_CSV(windowTimes, name_str, dayInterval):
    """Convert NumPy array of epoch times into human-readable format, in Pandas dataframe and CSV.
    
    Keyword argument:
    windowTimes -- NumPy array of dates in Unix epoch format (float)
    dayInterval -- Interval between windows in days (int)
    """
    windowDateTime = []
    dayInterval = str(dayInterval)
    
    for time in range(len(windowTimes)):
        # Convert ms to s by dividing by 1000
        windowDateTime.append( datetime.fromtimestamp(windowTimes[time] / 1000) )
        
    windowDateTime_df = pd.DataFrame(windowDateTime)
    windowDateTime_df.columns = ["Time"]
    windowDateTime_df.to_csv("windows_{}_{}day_interval.csv".format(name_str, dayInterval))

def ZonalStats(valueCollection, zoneCollection, statistic, scale = 10, tileScale = 16):
    """Computes zonal statistics across an image collection as a table. 
    An output value is computed:
    1) For every zone in the input zone collection.
    2) For every image in the value collection.
    
    Keyword arguments:
    valueCollection -- image collection whose values are used to find the output statistic
    zoneCollection -- feature collection of polygons that define zones to reduce the statistic to
    statistic -- Earth Engine reducer that calculates a given statistic for each zone
    scale -- pixel resolution of the data input into the function (default 10 m)
    tileScale -- sets operation tile size to avoid exceeding memory limit (default 16)
    """
    def ZS_Map(image):
        """Define zonal statistics operation, then apply to entire image collection
        
        Adapted from: https://gis.stackexchange.com/questions/333392/gee-reduceregions-for-an-image-collection
        """
        return image.reduceRegions(collection = zoneCollection,
                                   reducer = statistic,
                                   scale = scale,
                                   tileScale = tileScale)
    reduced = valueCollection.map(ZS_Map)
    
    # The above gives a collection of collections. To convert to a "table", simply flatten it.
    return reduced.flatten()

def MovingWindowDetectChange(pre_window_Start, window_Start, window_End, post_window_End, sizeWindows, numWindows, polygons, slope_threshold = 0.5, curv_threshold = -0.005):
    """Compute SAR intensity change detection algorithm over a moving window and aggregate with zonal statistics.
    The change detection uses the ratio of a stack of SAR images before the window and a stack of SAR images after the window.
    
    SAR intensity ratio change detection code written by Mong-Han Huang and Alexander L. Handwerger.
    https://doi.org/10.5194/nhess-2021-283
    
    Keyword arguments:
    pre_window_Start -- date to begin the pre-window stack
    window_Start -- date that the window begins
    window_End -- date that the window ends
    post_window_End -- date to end the post-window stack
    sizeWindows -- duration of the window in days; also determines moving window step size
    numWindows -- number of times to move the window forward
    polygons -- feature collection of polygons that define zones to detect landslides within
    slope_threshold -- upper threshold for slope mask (default 0.5 degrees)
    curv_threshold -- lower threshold for curvature mask (default -0.005 m/m^2)
    """
    import datetime as dt
    # Get change detection function, based on Handwerger et al. (2021)
    from ChangeDetect import I_Ratio
    
    # Define an area to filter image collections by.
    aoi = polygons.geometry().bounds()
    # Create an array to hold running list of window dates.
    windowEpochTime = np.empty(numWindows, dtype = float)
    
    ChangeCollection = []
    i = 1
    print("Processing . . .\n")
    
    # Run the change detection function a number of times while changing the window each iteration.
    for window in range(numWindows):
        ChangeImg = I_Ratio(aoi, slope_threshold, curv_threshold,
                            pre_window_Start, window_Start,
                            window_End, post_window_End)
        
        # Get the approximate date of the window by finding the mean of the beginning and end dates.
        windowDateAvg = (window_Start.getInfo()['value'] + window_End.getInfo()['value']) / 2
        
        # Divide by 1000 to convert from ms to s.
        windowDateStr = dt.datetime.fromtimestamp(windowDateAvg / 1000).strftime('%Y-%m-%d')
        print('\tWindow: ', i, ' of ', numWindows, '(',windowDateStr , ')')
        i = i + 1
        
        if ChangeImg is None:
            """Results from no data in one or more stacks.
            Pass null output and advance to the next window."""
            # Prepare for the next computation.
            # Move all dates forward a set number of days equal to the window size.
            pre_window_Start = pre_window_Start.advance(sizeWindows, 'day')
            window_Start = window_Start.advance(sizeWindows, 'day')
            window_End = window_End.advance(sizeWindows, 'day')
            post_window_End = post_window_End.advance(sizeWindows, 'day')
            
            windowEpochTime[window] = None
        else:
            # Add the date of this window to the list.
            windowEpochTime[window] = windowDateAvg
            
            pre_window_Start = pre_window_Start.advance(sizeWindows, 'day')
            window_Start = window_Start.advance(sizeWindows, 'day')
            window_End = window_End.advance(sizeWindows, 'day')
            post_window_End = post_window_End.advance(sizeWindows, 'day')

            # Build an image collection out of the intensity change images.
            # Saving all images to a list and converting all at once tends to run out
            # of memory, so add images to collection one at a time.
            if not ChangeCollection:
                # Initialize the collection during the first iteration of the loop.
                ChangeCollection = ee.ImageCollection(ChangeImg)
            else:
                # There is no EE method to add an image to an image collection, so
                # a dummy collection is needed to merge with the existing collection.
                ChangeCollection = ChangeCollection.merge(ee.ImageCollection(ChangeImg))
        
    # Find zonal statistics across the entire image collection.
    zonalChange_sum = ZonalStats(ChangeCollection, polygons, ee.Reducer.sum())
    zonalChange_count = ZonalStats(ChangeCollection, polygons, ee.Reducer.count())

    print('\n*** Complete! ***')
    return zonalChange_sum, zonalChange_count, windowEpochTime

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


# Name of output Google Drive folder for CSVs and GeoTiffs.
out_dir = 'EE_SAR_MovingWindow_ChangeDetection_OR_Eugene_Subset'
# Descriptive prefix for output files.
file_name = 'OR_LS_Eugene'

# Area of interest to make a grid over.
ee_import = ee.FeatureCollection("users/bvoelker/OR_Landslides/OR_LS_Eugene_Subset")
aoi = ee_import.geometry()

# The polygons to check for changes in will be a regularly spaced grid.
grid_size = 500 # 500m scale
grid = MakeGrid(aoi, grid_size) 
grid = Add_ID_to_Features(grid)
# Define duration of window.
sizeWindows = 7 # days

# Number of windows controls number of loops in change detection function.
# sizeWindows * numWindows = total detection period
numWindows = 208 + 20 # approx. 4 years plus 5 months

# Define how many days are in the pre window stack.
# Recommended to be multiples of one year to encompass whole cycles of seasonal vegetation change.
sizePreStack = 365 # days

# Define how many days are in the post window stack.
# Shorter windows will resolve changes more quickly at the cost of more noise in the SAR stack.
sizePostStack = 60 # days

# These parameters determine only the dates of the first window. All dates will be iteratively
    # moved forward within the 'MovingWindowDetectChange' function.
# The date to start initial pre-window stack.
"""GEE S1_GRD data beings on 2014-10-03T00:00:00Z.
First ASF S1 data over Oregon is June 3 2015.
However, full coverage of the state only begins in 2016."""
pre_window_Start = ee.Date('2016-06-03T00:00') # format: 'yyyy-mm-dd-HH':MM

# Date of the initial window start | End of pre-window stack.
window_Start = pre_window_Start.advance(sizePreStack, 'day')  # format: 'yyyy-mm-dd-HH':MM

# Date of the initial window end | Start of post-window stack.
window_End = window_Start.advance(sizeWindows, 'day')

# The date to end initial post-window stack.
post_window_End = window_End.advance(sizePostStack, 'day')

# Apply parameters to moving window SAR intensity change detection function.
sumChangePerZone, countPerZone, windowTimes = MovingWindowDetectChange(pre_window_Start, window_Start,
                                                                            window_End, post_window_End,
                                                                            sizeWindows, numWindows,
                                                                            grid)

# Remove nans from window date array, in cases where there were empty stacks.
windowTimes = windowTimes[~np.isnan(windowTimes)]

# To further manipulate the data outside of EE, export to a CSV.
exportToCSV_sum = ee.batch.Export.table.toDrive(collection = sumChangePerZone,
                                            folder = out_dir,
                                            description = "{}_99th_Ptile_sum".format(file_name),
                                            fileFormat = 'CSV')
exportToCSV_count = ee.batch.Export.table.toDrive(collection = countPerZone,
                                            folder = out_dir,
                                            description = "{}_pixel_count".format(file_name),
                                            fileFormat = 'CSV')
                                            
exportToCSV_sum.start()
exportToCSV_count.start()

WindowDate_to_CSV(windowTimes, file_name, sizeWindows)

"""
Wait for export to finish in GEE before continuing!
"""
#########################################################################################################################
"""
When finished, save CSVs to same directory as notebook.
"""

detection_date_list = []
detection_ID_list = []
detection_index_list = []

# Number of landslide polygons or grid cells.
numPolys = grid.size().getInfo()

# To avoid repead detections from the same event, determine an amount of time before
# another detection can be acquired.
reset_interval = math.ceil(math.ceil(sizePreStack / 2) / sizeWindows)

# Load in the CSVs:
# The needed columns are the ID numbers and zonal statistics result (either sum or count).
data_sum_raw = pd.read_csv("{}_99th_Ptile_sum.csv".format(file_name))
data_count_raw = pd.read_csv("{}_pixel_count.csv".format(file_name))

# Read in window dates
windows = pd.read_csv("windows_{}_{}day_interval.csv".format(file_name, sizeWindows))
windows['Time'] = pd.to_datetime(windows['Time'])

# Create an array of each polygon's unique identifier.
ids = data_sum_raw['ID'][0:numPolys]

# Format the mean intensity ratio:
# Rows: Zones (grid cells).
# Columns: Sum or count at each moving window date.
data_sum = np.transpose(data_sum_raw['sum']
                        .values
                        .reshape(-1, numPolys)
                       )
data_count = np.transpose(data_count_raw['count']
                          .values
                          .reshape(-1, numPolys)
                         )

# Dividing by a count of 0 results in NaN and a warning. Supress the warning.
np.seterr(all="ignore")
# Diving by the pixel count 'normalizes' data by the area.
# Practically it finds the percentage of pixels above 99th percentile.
data_normalized = data_sum / data_count
# Change NaN values to 0.
data_normalized = np.nan_to_num(data_normalized)


# Find the peaks in each time series.
for row_number, row in enumerate(data_normalized):
    """Adapted from: https://pythonawesome.com/overview-of-the-peaks-dectection-algorithms-available-in-python/"""
    peak_index = scipy.signal.argrelextrema(row,
                                            comparator = np.greater,
                                            order = 3)
    # Record possible detections according to three criteria:
    # 1) Peak is detected
    # 2) Above pixel percentage threshold of 10%
    # 3) More than half the time of the pre-window stack has passed.
    detect = False
    # Initialize as an arbitrary low number that is lower than first index.
    detect_index = -2*reset_interval
    for index in peak_index[0]:
        #Check if enough time has passed to reset the detection. If not, pass.
        if index < detect_index + reset_interval:
            pass
        elif row[index] > 0.1:
            """With a percentage threshold of 10%,
            Record the date and time of the detected peak."""
            detection_date_list.append(windows['Time'][index])
            detection_ID_list.append(int(ids[row_number]))
            detect = True
            detect_index = index

# Get list of indices too, which start at 0 (offset of 1 from the IDs)
#detection_index_list = [x - 1 for x in detection_ID_list]

# Preview the data.
#print(numPolys)
#print(len(detection_ID_list))
#print(detection_date_list)
#print(detection_ID_list)

# Depending on the data, one may need to change the size of the figure to suit the plot as needed.
figsize = (20,8)

"""From https://www.geeksforgeeks.org/python-ways-to-remove-duplicates-from-list/"""
unique_IDs = []
unique_indexes = []
[unique_IDs.append(x) for x in detection_ID_list if x not in unique_IDs]

# Get list of unique indices too, which start at 0 (offset of 1 from the IDs)
unique_indexes = [x - 1 for x in unique_IDs]

numDetections = len(unique_IDs)

# Create a square plot of subplots, based on the number of detections in the data.
plt_rows = math.ceil(math.sqrt(numDetections))
# Could use floor instead of ceil for the columns, but it doesn't work at very small plot sizes.
plt_cols = math.ceil(math.sqrt(numDetections))
fig, axes = plt.subplots(plt_rows, plt_cols, sharex=True, sharey=True, figsize = figsize)

for i, ax in enumerate(axes.flatten()):
    # There is a square array of subplots. But if there are not enough
    # plots to fill every cell, stop plotting when the limit is reached.
    if i > numDetections - 1:
        break
    
    currentData = data_normalized.T[:,unique_indexes[i]]
    newLabel = "Cell #" + str(unique_IDs[i])
    
    ax.plot(windows['Time'], currentData*100, color='black', label = newLabel)
    ax.axhline(0.05, alpha = 0.2, color = 'k', linestyle = '--')
    
    """Adapted from: https://pythonawesome.com/overview-of-the-peaks-dectection-algorithms-available-in-python/"""
    peak_index = scipy.signal.argrelextrema(currentData,
                                            comparator = np.greater,
                                            order = 3)
    ax.plot(windows['Time'][peak_index[0]], currentData[peak_index[0]]*100, 'ro')
    
    ax.set_xlabel('Date')
    ax.set_ylabel('% Pixels > 99th Ptile')
    ax.legend()
plt.savefig("{}_timeseries.png".format(file_name))

plt.show()

# Saves images to Google Drive.
[Get_BeforeAfter_Imagery(out_dir, ID, date, grid, sizeWindows) for ID, date in zip(detection_ID_list, detection_date_list)];

pre_NDVI = []
post_NDVI = []
for detection in range(len(detection_ID_list)):
    detID = detection_ID_list[detection]
    detDate = detection_date_list[detection]
    
    new_pre_NDVI, new_post_NDVI = Get_BeforeAfter_NDVI(detID, detDate, grid, sizeWindows)
    pre_NDVI.append(new_pre_NDVI.getInfo()['NDVI'])
    post_NDVI.append(new_post_NDVI.getInfo()['NDVI'])

NDVI_df = pd.DataFrame(list(zip(detection_ID_list, detection_date_list, pre_NDVI, post_NDVI)), columns = ['ID', 'Date', 'Pre NDVI', 'Post NDVI'])
NDVI_df.to_csv("{}_NDVI.csv".format(file_name))
print(NDVI_df)
