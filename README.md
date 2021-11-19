# Purpose
These scripts use Google Earth Engine to process SAR data from the ESA Sentinel-1 satellite platform for the purpose of change detection of the Earth's surface. Two stacks of SAR intensity images are created that cover two different periods of time. The time period between the stacks is effectively the 'detection window' where, if a change occurs, it will be detected by the algorithm. 

The median of each stack is found, which reduces noise and creates two images: one from the 'before' stack, and one from the 'after' stack. The ratio of the before and after stacks is then calculated. Since GEE S1 GRD data is in decibels, simply subtracting the two images is an equivalent operation to taking the ratio.

The ratio of the SAR intensity images is used to find areas of change -- larger changes in ground surface properties will create a larger ratio. This creates an intensity ratio image. Then a percentile threshold of the intensity ratio image is applied, where any pixels higher than the threshold are classified as a change. Handweger et al. (2021) found that for landslide detection, using the 99th percentile minimized false positives.

The input area of interest (Earth Engine geometry) is used to create a grid of rectangular polygon tiles. For each tile, the percentage of pixels over the value 99th percentile is calculated.

Next, this change detection process is repeated by moving the detection window and the dates of the stacks forward by a set amount of time. The result is a time series, for each cell, of percent area cover by pixels over the 99th percentile value.

Finally, for each time series, peak detection is applied to find points where the time series had a sharp increase, which corresponds to a time when the intensity ratio rises noticably. For efficiency, repeat detections from the same event are ignored and a baseline percent pixel threshold is applied to avoid false detections from random noise. A list of all cells with detections, and the date of detections, is recorded. This list is used to generate an optical satellite image of the area before/after the event, and the average NDVI of the area before/after the event. 

# How to use

1) Run Iratio_ChangeDetect_GenCSV.py
2) Wait for Google Earth Engine to export the two CSVs at the end of the code (could take more than 1 hour).
3) Run Iratio_ChangeDetect_Process_CSV.py

Requires use of an Earth Engine account.

If GEE export fails, it is likely because there are no images at the beginning of the time series. Try a later start date.

Recommended polygon tile spacing is 500m. Earth Engine may run out of memory with smaller tiles.

Variables may be changed to suit desired area of interest. For example, input different Earth Engine geometry/asset to run over a different area. Or change the intitial time and step sizes to run the time series over different intervals. Ensure that any variables shared between the two Python scripts are the same.

***********
ChangeDetect.py  contains a SAR intensity ratio change detection function which is adapted from Handweger et al. (2021):

Handwerger, A. L., Jones, S. Y., Amatya, P., Kerner, H. R., Kirschbaum, D. B., and Huang, M.-H.: Strategies for landslide detection using open-access synthetic aperture radar backscatter change in Google Earth Engine, Nat. Hazards Earth Syst. Sci. Discuss. [preprint], https://doi.org/10.5194/nhess-2021-283, in review, 2021.

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
