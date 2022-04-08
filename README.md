InLineAnalysis
==============

InLineAnalysis provide a set of routine to process observations from scientific instrumentation collected under the way of research vessels. The code can be used to synchronize the stream of data from each instrument; separate the different periods of water going through the instruments (Filtered, Total, DIW); average the data every minute; and run manual and automatic quality check (QC); calibrate, adjust, and compute products from the raw data. The output can then be exported to SeaBASS format to share it with the community.

To date the following instruments are supported:
  + TSG
  + ACS
  + AC9
  + FLBBCD
  + BB3
  + WSCD
  + LISST
  + ALFA
  + PAR
  + SUVF
  + HyperBB
  + NMEA (GPS)
  + Atlas TSG

The application is under continuous development, please use the [GitHub Issues Tracker](https://github.com/OceanOptics/InLineAnalysis/issues) to signal any bug, or feature requests. The documentation is limited to this files and the comments in the code at this time.

# Quick start
## Installation
Set the matlab working folder to the InLineAnalysis folder:
  
    cd ~/Documents/MATLAB/InLineAnalysis

Create the packages folder
  
    mkdir packages

Download and add the packages below in the `packages` folder:
  + [Datetick Zoom, Pan & Subplot with Day of Year](https://www.mathworks.com/matlabcentral/fileexchange/25927-datetick-zoom--pan---subplot-with-day-of-year)
  + [Spectral and XYZ Color Functions](https://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions?focused=5172034&tab=function)
  + [Progress bar](https://github.com/JAAdrian/MatlabProgressBar)
  + [Scattering in pure seawater betasw_ZHH2009.m](https://github.com/ooici/ion-functions/blob/master/ion_functions/data/matlab_scripts/flort/betasw_ZHH2009.m)
  + [TOES-10 Toolbox](https://www.teos-10.org/software.htm)

Install the [JSONlab toolbox](https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files)

Recommended Toolbox
  + image processing toolbox
  + Finance
  + Parallel

## Configuration
For each project a new configuration file should be written. The configuration files are located in the `cfg` directory.

## Run step by step
Adapt the `main.m` script to process your InLine data. To consolidate the dataset and export it, adapt the script `Consolidate.m`.

# Code description
The application is accessed using the `InLineAnalysis` class.
The `Instrument` class provide a common interface for each instrument. 
The `lib` folder contains the core of the application, most of those function can be used separately.
