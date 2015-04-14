#!/bin/bash
. /reg/g/psdm/etc/ana_env.sh
. /reg/g/psdm/bin/sit_setup.sh /reg/neh/home/alindahl/amoi0114/analysis/ana-0.13.5

# PV rellated libraries are found here
export PATH=$PATH:/reg/common/package/release/sxr-0.0.1/x86_64-rhel5-gcc41-opt/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/reg/common/package/release/sxr-0.0.1/x86_64-rhel5-gcc41-opt/lib
