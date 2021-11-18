# How to use:

1) Run Iratio_ChangeDetect_GenCSV.py
2) Run Iratio_ChangeDetect_Process_CSV.py

Requires use of an Earth Engine account.

If GEE export fails, it is likely because there are no images at the beginning of the time series. Try increasing the

Variables may be changed to suit desired area of interest. For example, input different Earth Engine geometry/asset to run over a different area. Or change the intitial time and step sizes to run the time series over different intervals.

***********
ChangeDetect.py  contains a function which is adapted from Handweger et al. (2021):

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
