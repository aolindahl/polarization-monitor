# polarization-monitor
Polarization monitor using cookie box

# Setup
polarization-monitor includes aolPyModules as a git submodule. In order to correctly initialize the submodule run the following commands in the polarization-monitor directory.
```
git submodule init
git submodule update
```

# Description
The main script is the onlineDataCruncher.py which can be run on a single core or on multiple cores using mpirun. Refer to the built in help (python onlineDataCruncher.py --help) for available command line options.

onlineDataCruncher.py uses cookie_box.py from aolPyModules which in turn relies on tof.py (also in aolPyModules). A configuration file (rally a python module) is loaded and passed to the CookieBox object. The configuration file defines many of the properties of the analysis. The goal is that all options and settings should be available at the command line or in the config file. The config file to be used should be specified using the -c option.
In the config file the settings for each of the 16 TOF detectors in the cooke box are specified. One of these settings defines the time to energy conversion as a path to a json file containing a definition of conversion parameters.

# Running instructions
login as amoopr on the right daq-amo mon machines (see amo.cnf)
```
cd /reg/neh/operator/amoopr/experiments/polarization-monitor
```
get psana environment:
```
source /reg/g/psdm/etc/ana_env.sh
cd ~amoopr/experiments/polarization-monitor
```
command to run on a mon node (e.g. daq-amo-mon02) to see
what names of devices are in the daq:
```
psana -m EventKeys -n 1 shmem=AMO.0
```
run on 1 core:
```
python onlineDataCruncher.py -v -c config15-04-07.py -tA 100 -1A 10
```
run in parallel on many 3 nodes (24 cores):
```
`which mpirun` -n 24 --host daq-amo-mon02,daq-amo-mon03,daq-amo-mon04 ./onlineExecute.sh
```
run plots (e.g. on amo-console machine as user amoopr):
```
python plotter.py -H daq-amo-mon02
```
