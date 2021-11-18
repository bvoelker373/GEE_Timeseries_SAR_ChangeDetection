# Initialize Earth Engine Python API.
import ee
ee.Authenticate()
ee.Initialize()

# For data manipulation and analysis.
import pandas as pd
import numpy as np
np.set_printoptions(precision=4, suppress=True)
from datetime import datetime

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

def WindowDate_to_DFandCSV(windowTimes, dayInterval):
    """Convert NumPy array of epoch times into human-readable format, in Pandas dataframe and CSV."""
    windowDateTime = []
    dayInterval = str(dayInterval)
    
    for time in range(len(windowTimes)):
        # Convert ms to s by dividing by 1000
        windowDateTime.append( datetime.fromtimestamp(windowTimes[time] / 1000) )
        
    windowDateTime_df = pd.DataFrame(windowDateTime)
    windowDateTime_df.columns = ["Time"]
    windowDateTime_df.to_csv("windows_{}day.csv".format(dayInterval))
    return windowDateTime_df

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
    # Get change detection function, based on Handwerger et al. (2021)
    from ChangeDetect import I_Ratio
    
    # Define an area to filter image collections by.
    aoi = polygons.geometry().bounds()
    # Create an array to hold running list of window dates.
    windowEpochTime = np.empty(numWindows, dtype = float)
    
    # Run the change detection function a number of times while changing the window each iteration.
    for window in range(numWindows):
        ChangeImg = I_Ratio(aoi, slope_threshold, curv_threshold,
                            pre_window_Start, window_Start,
                            window_End, post_window_End)
        
        # Get the approximate date of the window by finding the mean of the beginning and end dates.
        windowDateAvg = (window_Start.getInfo()['value'] + window_End.getInfo()['value']) / 2
        windowEpochTime[window]= windowDateAvg
        # Prepare for the next computation.
        # Move all dates forward a set number of days equal to the window size.
        pre_window_Start = pre_window_Start.advance(sizeWindows, 'day')
        window_Start = window_Start.advance(sizeWindows, 'day')
        window_End = window_End.advance(sizeWindows, 'day')
        post_window_End = post_window_End.advance(sizeWindows, 'day')
        
        # Build an image collection out of the intensity change images.
        # Saving all images to a list and converting all at once tends to run out
        # of memory, so add images to collection one at a time.
        if window == 0:
            # Initialize the collection during the first iteration of the loop.
            ChangeCollection = ee.ImageCollection(ChangeImg)
        else:
            # There is no EE method to add an image to an image collection, so
            # a dummy collection is needed to merge with the existing collection.
            ChangeCollection = ChangeCollection.merge(ee.ImageCollection(ChangeImg))
            
        
    # Find zonal statistics across the entire image collection.
    zonalChange_sum = ZonalStats(ChangeCollection, polygons, ee.Reducer.sum())
    zonalChange_count = ZonalStats(ChangeCollection, polygons, ee.Reducer.count())
    
    return zonalChange_sum, zonalChange_count, windowEpochTime


# Output Drive folder for CSVs and GeoTiffs.
out_dir = 'EE_SAR_MovingWindow_ChangeDetection_HolyFire'

# Area of interest to make a grid over.
ee_import = ee.FeatureCollection("users/bvoelker/Wildfires/Holy_Fire_2018_08_06_Boundary")
#ee_import = ee.FeatureCollection("users/bvoelker/hiroshima_landslide_subset2")
aoi_bounding_box = ee_import.geometry().bounds();

# The polygons to check for changes in will be a regularly spaced grid.
grid_size = 500 # 500m scale
grid = MakeGrid(aoi_bounding_box, grid_size) 
grid = Add_ID_to_Features(grid)

# Define duration of window.
sizeWindows = 7 # days
# Number of windows controls number of loops in change detection function.
# sizeWindows * numWindows = total detection period
numWindows = 52
# Define how many days are in the pre window stack.
# Recommended to be multiples of one year to encompass whole cycles of seasonal vegetation change. (?)
sizePreStack = 365 # days
# Define how many days are in the post window stack.
# Recommended to be shorter, e.g. 1-3 months to capture change rapidly with few 'repeat' detections. (?)
sizePostStack = 60 # days

# The date to start initial pre-window stack.
pre_window_Start = ee.Date('2018-08-06T00:00') # format: 'yyyy-mm-dd-HH':MM
# Date of the initial window start.
window_Start = pre_window_Start.advance(sizePreStack, 'day')  # format: 'yyyy-mm-dd-HH':MM
# Date of the initial window end.
window_End = window_Start.advance(sizeWindows, 'day')
# The date to end initial post-window stack.
post_window_End = window_End.advance(sizePostStack, 'day')


# Apply parameters to moving window SAR intensity change detection function.
sumChangePerZone, countPerZone, windowTimes_7day = MovingWindowDetectChange(pre_window_Start, window_Start,
                                                                            window_End, post_window_End,
                                                                            sizeWindows, numWindows,
                                                                            grid)
# To further manipulate the data outside of EE, export to a CSV.
exportToCSV_sum = ee.batch.Export.table.toDrive(collection = sumChangePerZone,
                                            folder = out_dir,
                                            description = 'HF_CA_99th_Ptile_sum',
                                            fileFormat = 'CSV')
exportToCSV_count = ee.batch.Export.table.toDrive(collection = countPerZone,
                                            folder = out_dir,
                                            description = 'HF_CA_pixel_count',
                                            fileFormat = 'CSV')


"""Start Earth Engine export."""
exportToCSV_sum.start()
exportToCSV_count.start()

# Export the window dates to a CSV.
window_df = WindowDate_to_DFandCSV(windowTimes_7day, 7)
