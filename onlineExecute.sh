#!/bin/bash
cd ~amoopr/experiments/polarization-monitor

source /reg/g/psdm/etc/ana_env.sh
source /reg/g/psdm/bin/sit_setup.sh
source pv_environment_setup.sh

python onlineDataCruncher.py -v -c config15-04-07.py -tA 100 -1A 10
