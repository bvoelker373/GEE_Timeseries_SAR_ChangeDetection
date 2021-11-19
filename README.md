# Purpose
Ratio of two SAR backscattered intensity images can be used to detect changes of ground surface properties.
Create stacks of many SAR images before and after event to reduce noise.
Moving the time period between images (detection window) allows time series to be constructed.
Finding anomalies or peaks in time series could pinpoint landslide events of interest.
Use Google Earth Engine to handle large volumes of data.
![image](https://user-images.githubusercontent.com/94650022/142574427-22312724-5cc3-43f1-8fb9-e730d5e1e2ec.png)

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
