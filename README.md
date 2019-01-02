the basic analysis framework for t576.

TODO:  fix the positions for some remaining events (run 4 mostly). 

# Dependencies:

The only dependency for this software is ROOT:

https://root.cern.ch/downloading-root

and i recommend installing the latest version. (6.14.06 at time of writing). installing root6 is very straightforward thanks to CMAKE.

 #### be sure, once installed, that you do the "source /path/to/root/bin/thisroot.sh" so that your environment variables are set correctly.



additionally, this software uses the very nice, lightweight code for reading python .npz and .npy files into c++ ```cnpy```, which can be found at

https://github.com/rogersce/cnpy

we have used the relevant source from cnpy in this project, and lisence info can be found in LICENSE. 

# Install:

the software uses CMAKE, and you need to set a couple environment variables. T576_INSTALL_DIR is the top directory for the install. headers will install to T576_INSTALL_DIR/include/t576, sources to  T576_INSTALL_DIR/src/t576,  run logs and index files to T576_INSTALL_DIR/share/t576 and the t576 library to T576_INSTALL_DIR/lib. any directories below T576_INSTALL_DIR/ which don't exist will be automatically created.

T576_DATA_DIR is where the data is, more on that below.

to install: 

first, add these to your bashrc:
```bash

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


then, (optional) tell it where you'd like to install the software, as above. if you don't do this, it will install to /usr/local/.

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

double z = ev->scope->pos[2].Z();

//these event-by-event variables are all loaded up when you call
//ev->loadScopeEvent(event number). 

//the same can be done for the surf, of course
ev->surf->loadSurfEvent(9000);

//which will put all the things in the proper place, e.g. to draw a graph
//for channel 9 in the surf for this event, do

ev->surf->gr[9]->Draw("al");

//furthermore, you can set the level of interpolation for all events

ev->setInterpGsNs(20);

//which will set the sample rate for all resultant graphs to 20GS/s.
//subsequent calls to ev->loadSurf[Scope]Event(xx) will have this level
//of interpolation applied. details are in the source.


```
future releases will have all sorts of things, like pointing maps and correlations etc.

# Data:

this will look for data in the directory you specify with T576_DATA_DIR. if you have a full copy of the data, it should be in the format

/path/to/t576/run2/(root/ py/) where root/ and py/ are under run2/. therefore, for this structure, you'd put

```bash
export T576_DATA_DIR=/path/to/t576/run2
```
in your .bashrc. this will let the program find the data.

#### filenames:

all of the filenames are stored with a timestamp then the run major and the run minor, e.g. 20181031122345run1_999.root where 1 is the major, 999 is the minor. the timestamp just helps to sort the files. the major and minor are used to index in root. 

#### missing data:

there is a data flag, isGood(), which says whether an event is good. it is good if all of the event info is present, and if there wasn't a note in the run log for the event being corrupted in some way.

if data was missing, like we don't have a distance in place or something, the value is 999. so if event->surf{scope}->pos.Mag()=999, the position is not known. 

# Problems:

are going to happen. email me.


have fun!