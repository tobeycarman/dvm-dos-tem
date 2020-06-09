#!/bin/bash

# setting up to use boost on modex
echo "Loading module files..."
echo "NOTE: This will work with netcdf/4.3.3.1 or netcdf/4.7.1-gnu540" 
module purge; module load python/3.6.2 gcc/5.4.0 boost/1.67.0 netcdf/4.7.1-gnu540 jsoncpp/jsoncpp-1.8.4

echo "Setting up site specific inlcudes..."
export SITE_SPECIFIC_INCLUDES="-I/data/software/src/jsoncpp_1.8.4/jsoncpp-1.8.4/include/ -I/data/software/src/openblas/OpenBLAS-0.3.7/lapack-netlib/LAPACKE/include/"

echo "Setting up site specific link flags..."
echo "NOTE: I have no idea why this has become an issue all of a sudden. Seems "
echo "      to have to do with linking openblas, but for some reason now we get "
echo "      a netcdf error about unresolved symbols in TEMUtilityFunctions.cpp?!?"
echo "      Interwebs say that this can be an issue if you don't have netcdf 4 "
echo "      enabled but I checked that we do (nc-config). And this hasn't been an "
echo "      issue on any other machine with dvmdostem v0.2.3. So here we add some "
echo "      special options to the gcc link step."
echo ""
export SITE_SPECIFIC_LINK_FLAGS="-Wl,--unresolved-symbols=ignore-in-object-files"

echo "Setting up site specific libs..."
export SITE_SPECIFIC_LIBS="-L/data/software/boost/1.67.0/lib/ -L/data/software/src/jsoncpp_1.8.4/jsoncpp-1.8.4/build-shared/ -L/data/software/src/openblas/OpenBLAS-0.3.7/"

echo "Adjusting Makefile to use openblas (for lapacke libary)."
sed -e 's:-llapacke:-lopenblas:' Makefile > Makefile.tmp && mv Makefile.tmp Makefile

echo "Setting LD_LIBRARY_PATH to pick up openblas at run time."
export LD_LIBRARY_PATH="/data/software/src/openblas/OpenBLAS-0.3.7/:$LD_LIBRARY_PATH"

echo "NOTE: This file will NOT work if it is run as a script!"
echo "      Instead use the 'source' command like this:"
echo "      $ source env-setup-scripts/setup-env-for-atlas.sh"
echo ""
echo "NOTE: Please remember not to commit the modified Makefile!"
echo "      You can revert the change with this command:"
echo "      $ git checkout -- Makefile"


