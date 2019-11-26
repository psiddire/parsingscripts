# parsingscripts

## Loading fastjet 

The fastjet directory in this repository is already compiled. It can be loaded in the python scripts by adding it's path to bash_profile:

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'pwd'/fastjet-install/lib
PYTHONPATH=$PYTHONPATH:'pwd'/fastjet-install/lib/python3.7/site-packages

In the above lines you need to edit the 'pwd' to the directory path.

In case there is some error with the above procedure you can re-build the fastjet again using:

./configure --prefix='pwd'/../fastjet-install --enable-pyext \s\s
make \s\s
make check \s\s
make install \s\s


## Background

For information on the LHC Olympics 2020, see [this website](https://indico.cern.ch/event/809820/page/16782-lhcolympics2020).


