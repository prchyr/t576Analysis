the basic analysis framework for t576.

TODO:  fix the positions for some remaining events (run 4 mostly). expand the TUtil class to include the SVD.


# Dependencies:

there are 2 dependencies for this software:

1) ROOT:

https://root.cern.ch/downloading-root

and i recommend installing the latest version. (6.14.06 at time of writing). installing root6 is very straightforward thanks to CMAKE.

 #### be sure, once installed, that you do the "source /path/to/root/bin/thisroot.sh" so that your environment variables are set correctly.

2) CLHEP:

clhep is a really great library for doing physics things with code. it has datatypes like 3-vectors, 4-vectors, and a bunch of physical constants and a nice system of units that makes coordinate transformations and relativistic things like boosts a breeze.

https://proj-clhep.web.cern.ch/proj-clhep/clhep23.html

and again i recommend installing the latest version (2.4.1 at time of writing). installation is pretty straightforward for this too.


# Install:

the software uses CMAKE, and you need to set a couple environment variables.T576_INSTALL_DIR is the top directory for the install, under which must be /include /src /share and /lib. (if those aren't there, they will be made.) T576_DATA_DIR is where the data is, more on that below. 2 example installs follow.

to install: 

first, add these to your bashrc:
```bash
#point to the build directory for CLHEP, which conatains CLHEPConfig.cmake
export CLHEP_DIR=/path/to/CLHEP/build/dir
#point to where you want the tools installed. default is /usr/local
export T576_INSTALL_DIR=/path/to/your/fav/install/dir
#point to the t576 data. 
export T576_DATA_DIR=/path/to/t576/data
```

then,
```bash
cd /your/fav/source/directory
git clone https://github.com/prchyr/t576Analysis.git
cd t576Analysis

./install.sh
```

and it will compile, install, and then compile a little test macro that you can run to see that everything worked correctly.


details about the above:


Firstly for anything to work you need to tell it where to find ROOT. this is usually done already when you run the "source thisroot.sh" script from root. I'll assume you know what this is already. to know if things are sourced correctly, type "echo $ROOTSYS", and you should get the correct directory to ROOT.

next, tell it where CMAKE can find CLHEP. you can do that in 2 ways. you can either add the installation path for CLHEP to your CMAKE_PREFIX_PATH, or you can point CMAKE to the CLHEP build directory. this is the build directory you made when you installed CLHEP. it can be inside of the CLHEP directory, or somewhere else.


then, (optional) tell it where you'd like to install the software. if you don't do this, it will install to /usr/local/, meaning that it will put the header files in /usr/local/include, the t576 library into /usr/local/lib, and the share files into /usr/local/share. If you set your own install directory, /your/fav/install/dir, headers will go to /your/fav/install/dir/include, the library will go to /your/fav/install/dir/lib, etc. 

finally, tell it where the t576 data lives. the data directory structure must be /path/to/t576/data/(root/ py/) meaning root/ and py/ live under some common directory. and under py/, lives dat/ and ped/. contact me if there are questions on this.

# Using:

to use inside of root is simple, just point ROOT to the t576 library that you just made on installation. i recommend doing this in your rootlogon.C file by adding this line like so:
```c++
gROOT->ProcessLine(".L $T576_INSTALL_DIR/lib/libt576.so");
```
or if you didn't set $T576_INSTALL_DIR, do
``` c++
gROOT->ProcessLine(".L /usr/local/lib/libt576.so");
```


to use the software in your own standalone programs, just include the header and compile against the root libraries and the t576 library like so:

```c++
#include "t576/T576Event.hh"
```
and compile with
```bash
g++ -std=c++11 your_script.cc -o your_script `root-config --cflags --glibs --libs` -L$T576_INSTALL_DIR/include/t576 -lt576
```
ignore the -L flag if you have $T576_INSTALL_DIR/include in your include path.


# T576Event class

this is the main workhorse. you can load an event by run number or by an overall index number, which is generated the first time you use the class. if you run the example that is built on installation, this will be built at that time.

```c++
//make a new T576 event

T576Event * ev = new T576Event();

//load an event by run major, run minor, and event number within that file
//NOTE: you don't need to have the exact filename. they have timestamps etc
//in the names, but this will find it just based on the major and minor.

ev->loadScopeEvent(4, 6, 1);

ev->scope->gr[3]->Draw("al");

//will draw channel 4 on ths scope (ICT) as a graph.

//then you can try loading an event a different way:
ev->loadScopeEvent(27000);

//which will load whatever was the 27,000th event taken by the scope.
//to find out what run that was in,

cout<<ev->major<<" "<<ev->minor<<endl;
//or
cout<<ev->scopeFilename->Data()<<endl;

//there are other useful things like
ev->charge;
//which gives you the charge from the ICT. also
ev->scope->pos[2];
//would give you an Hep3Vector of the position of the channel 3 antenna,
//such that you can acces the z dimension like

double z = ev->scope->pos[2].z();

//these event-by-event variables are all loaded up when you call
//ev->loadScopeEvent(event number). 
```

# Data:

this will look for data in the directory you specify with T576_DATA_DIR. if you have a full copy of the data, it should be in the format

/path/to/t576/run2/(root/ py/) where root/ and py/ are under run2/. therefore, for this structure, you'd put

```bash
export T576_DATA_DIR=/path/to/t576/run2
```
in your .bashrc. this will let the program find the data.

FILENAMES: all of the filenames are stored with a timestamp then the run major and the run minor, e.g. 20181031122345run1_999.root where 1 is the major, 999 is the minor. the timestamp just helps to sort the files. the major and minor are used to index in root. 

# Problems:

are going to happen. email me.


have fun!