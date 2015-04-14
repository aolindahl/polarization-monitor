# polarization-monitor
Polarization monitor using cookie box

# Setup
polarization-monitor includes aolPyModules as a git submodule. In order to correctly initialize the submodule run the following commands in the polarization-monitor directory.

git submodule init
git submodule update

# Description
The main script is the onlineDataCruncher.py which can be run on a single core or on multiple cores using mpirun. Refer to the built in help (python onlineDataCruncher.py --help) for available command line options.

onlineDataCruncher.py uses cookie_box.py from aolPyModules which in turn relies on tof.py (also in aolPyModules). A configuration file (rally a python module) is loaded and passed to the CookieBox object. The configuration file defines many of the properties of the analysis. The goal is that all options and settings should be available at the command line or in the config file. The config file to be used should be specified using the -c option.
In the config file the settings for each of the 16 TOF detectors in the cooke box are specified. One of these settings defines the time to energy conversion as a path to a json file containing a definition of conversion parameters.
