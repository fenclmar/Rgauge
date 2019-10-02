# Rgauge
Functions for tipping bucket rain gauge calibration and data processing. Designed specifically for rain gauges and weather stations manufactured by FIEDLER AMS s.r.o. (https://www.fiedler.company/)

# Instalation
The repo is currently not ready to be directly compiled to R package. Install function by sourcing the lib_RGcalib.r file.

```
source(url('https://raw.githubusercontent.com/fenclmar/Rgauge/master/lib_RGcalib.r'))
```

The functions include comments which mostly provide a basic info what a function does and an info about the arguments and outputs
of a function. The functions are intended for these tasks:
- evaluation of laboratory (dynamic) calibraiton of rain gauges
- dynamic calibration of rain gauge data
- reading of logged data (in different formates)
- data cleaning and processing
- analysis of rain gauge data
