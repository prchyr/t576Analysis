#!/bin/bash

top_dir=$(pwd)

if [ -d "build" ]
then
    echo "deleting existing build directory, starting from scratch (it's ok)"
    rm -R "build"
fi

export T576_SOURCE_DIR=$top_dir

if [ -z "$1" ]
then
    if [ -z "$T576_INSTALL_DIR" ]
    then
	echo "No install directory specified. will install to /usr/local"
	export T576_INSTALL_DIR=/usr/local
    else
	echo "will install to $T576_INSTALL_DIR"
    fi
else
    export T576_INSTALL_DIR="$1"
    echo "will install to $T576_INSTALL_DIR"
    
fi



mkdir -p build &&
    cd build

cmake $top_dir &&
    make -B -j4 &&
    make install
