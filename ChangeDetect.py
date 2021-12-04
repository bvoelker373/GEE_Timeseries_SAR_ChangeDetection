"""
Converted from JavaScript code written by Mong-Han Huang and Alexander L. Handwerger.
https://doi.org/10.5194/nhess-2021-283

***********
MIT License

Copyright (c) 2021 Mong-Han Huang and Alexander L. Handwerger

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
***********
"""

def I_Ratio(aoi, slope_threshold, curv_threshold,
            preEventTime_1, preEventTime_2,
            postEventTime_1, postEventTime_2):
    """Calculates amplitude ratio of SAR stacks, adapted from Handwerger et al. (2020).
    Builds a stack of images before a landslide event and a stack after.
    The event must occur in the time window between the stacks.
    
    Keyword arguments:
    aoi -- geometry used to filter image collections to area of interest
    slope_threshold -- upper threshold for slope mask (degrees)
    curv_threshold -- lower threshold for curvature mask (m/m^2)
    preEventTime_1 -- start date of pre-event SAR stack 
    preEventTime_2 -- end date of pre-event SAR stack 
    postEventTime_1 -- start date of post-event SAR stack 
    postEventTime_2 -- end date of post-event SAR stack
    """
    import ee
    
    def LowEdgeMask(image):
        """Remove low edge values as suggested by GEE."""
        edge = image.lt(-30.0)
        maskedImage = image.mask().And(edge.Not())
        return image.updateMask(maskedImage)
    
    
    # Load Shuttle Radar Topography Mission (SRTM) Digital Elevation Model (DEM).
    dataset = ee.Image('USGS/SRTMGL1_003')
    elevation = dataset.select('elevation')
    slope = ee.Terrain.slope(elevation); #slope in degrees
    mask_slope = slope.gte(slope_threshold); # slope mask with values 0 or 1
    
    # Create water mask.
    waterZones = ee.Image('NASA/NASADEM_HGT/001').select('swb').clip(aoi)
    waterMask = waterZones.eq(0)
    
    # Define a Gaussian kernel to reduce noise in S1 scenes.
    # Uncomment here and below to use
    #smooth_S1 = ee.Kernel.gaussian(radius = 50, sigma = 20, units = 'meters', normalize = True)
    
    # Load Sentinel-1 (S1) amplitude data in VH polarization.
    imgVH = (ee.ImageCollection('COPERNICUS/S1_GRD')
                 .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                 .filter(ee.Filter.eq('instrumentMode', 'IW'))
                 .select('VH')
                 .filterBounds(aoi)
                 .map(LowEdgeMask))
    desc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')); #descending acquisition geometry data
    asc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));  #ascending acquisition geometry data
    
    # Calculate curvature.
    # Define a Gaussian kernel for smoothing. This step helps reduce noise in the curvature maps.
    smooth_curv = ee.Kernel.gaussian(radius = 120, sigma = 60, units= 'meters', normalize = True,)
    xyDemGrad = elevation.convolve(smooth_curv).gradient()
    xGradient = xyDemGrad.select('x').gradient()
    yGradient = xyDemGrad.select('y').gradient()
    curvature = xGradient.select('x').add(yGradient.select('y'))
    mask_curvature = curvature.gte(curv_threshold)
    
    
    #############################################################################
    ###                   Sentinel-1 SAR Image Processing                     ###
    #############################################################################
    preEventPeriod = ee.Filter.date(preEventTime_1,preEventTime_2)
    postEventPeriod = ee.Filter.date(postEventTime_1,postEventTime_2)
    
    # Median pre-event S1 SAR amplitude within the specified date ranges.
    preEventPeriod_asc = ee.Image.cat(asc
                                      .filter(preEventPeriod)
                                      .median()
                                     )
    preEventPeriod_desc = ee.Image.cat(desc
                                      .filter(preEventPeriod)
                                      .median()
                                      )
    
    # Post-event amplitude without filter.
    postEventPeriod_asc = ee.Image.cat(asc
                                       .filter(postEventPeriod)
                                       .median())
    postEventPeriod_desc = ee.Image.cat(desc
                                       .filter(postEventPeriod)
                                       .median())

    # Errors will arise from lack of SAR data for a given time period, and the stack will be empty.
    # Catch if any one of the stacks of empty and outpull a null value as the result.
    # The test is whether there are even any bands within the image created from the stacks.
    # Otherwise, if all the data is accounted for, proceed.
    if not all([preEventPeriod_asc.bandNames().getInfo(),
                postEventPeriod_asc.bandNames().getInfo(),
                preEventPeriod_desc.bandNames().getInfo(),
                postEventPeriod_desc.bandNames().getInfo()]):
        pass
    else:
        # Calculate the log ratio (using subtraction since data are in log scale)
        # for Pre- and Post-event S1 SAR Backscatter.
        Iratio_desc = preEventPeriod_desc.subtract(postEventPeriod_desc)
        Iratio_desc = Iratio_desc.clip(aoi)

        Iratio_asc = preEventPeriod_asc.subtract(postEventPeriod_asc)
        Iratio_asc = Iratio_asc.clip(aoi)

        Iratio_mean_desc_asc = (Iratio_asc.add(Iratio_desc)).divide(2)


        I_ratio_mean_masked = Iratio_mean_desc_asc.updateMask(mask_slope).updateMask(mask_curvature).updateMask(waterMask) 
        #I_ratio_mean_masked = Iratio_mean_desc_asc.updateMask(mask_slope).updateMask(mask_curvature)

        # Calculate percentiles of I_ratio. It is recommended to only use the 99th percentile.
        I_ratio_Percentiles = I_ratio_mean_masked.reduceRegion(
            reducer = ee.Reducer.percentile([90, 95, 99]),
            geometry = aoi,
            scale = 10,
            bestEffort = True
        )
        I_ratio_99thPtile = I_ratio_mean_masked.gte(ee.Number(I_ratio_Percentiles.get("VH_p99")))
        
        return I_ratio_99thPtile
