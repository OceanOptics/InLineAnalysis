InLineAnalysis
==============

This set of Matlab code is designed to process In-Line measurements done continuously on the flowthrough system of research vessel. It will synchronize the stream of data from each instrument; separate the different periods of water going through the instruments (Filtered, Total, DIW); bin the data ; and run an automatic quality check (QC) on the data as well as provide a graphical interface to manually QC the data. The processed data is saved in matlab format and can then be exported to SeaBASS format with minimal effort from the user.

To date the following instruments are supported:
  + TSG
  + ACS
  + FLBBCD
  + BB3
  + WSCD
  + LISST

The application is under continuous development and still need testing, please use the [GitHub Issues Tracker](https://github.com/OceanOptics/InLineAnalysis/issues) to signal any bug, or feature requests.
The documentation is minimal if not existing at the time, however, some help can be found in the comments of the code.

# Quick start
## Installation
Set the matlab working folder to the InLineAnalysis folder:
  
    cd ~/Documents/MATLAB/InLineAnalysis

Create the packages folder
  
    mkdir packages

Download and add the packages below in the `packages` folder:
  + [Datetick Zoom, Pan & Subplot with Day of Year](https://www.mathworks.com/matlabcentral/fileexchange/25927-datetick-zoom--pan---subplot-with-day-of-year)
  + [Spectral and XYZ Color Functions](https://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions?focused=5172034&tab=function)

Install the [JSONlab toolbox](https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files)

## Configuration
For each project processed a new configuration file must be written. The configuration files are located in the directory `cfg`, some example of configurations are available there.

## Run step by step
For each project the code that will be run should be written in `main_<project_name>.m`. Examples of the main script are available at the root of the project.

# Code description
The core of the application is in the `InLineAnalysis` class.
The `Instrument` class provide a common interface for each instrument. 
The `lib` folder contains the logic of the application, most of those function can be used separately.
