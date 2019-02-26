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
    cp $T576_SOURCE_DIR/example/loadEvents.cc build
    cd build


cmake $top_dir &&
    make -B -j4 &&
    make install &&
    printf "\n\n\n\na test script is being compiled...\n\n" &&

    g++ -std=c++11  loadEvents.cc -o loadEvents -I$T576_INSTALL_DIR/include `root-config --cflags --glibs --libs` -lt576 &&

    mv loadEvents $T576_INSTALL_DIR/share/t576/example/ &&
    printf "done. make sure that T576_INSTALL_DIR/include is in your include path, then try running:\n\n" &&
    printf " $T576_INSTALL_DIR/share/t576/example/loadEvents\n\n" &&
    printf "other examples have been built and now live in $T576_INSTALL_DIR/share/t576/example/ as well.\n" 


