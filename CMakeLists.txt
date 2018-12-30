# $Id: CMakeLists.txt 100923 2016-11-03 10:50:34Z gcosmo $

#



#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(t576)



unset(CMAKE_MODULE_PATH)
#use the CLHEP things
#will generate a bunch of library modules CLHEP::stuff that can be linked to
find_package(CLHEP REQUIRED)
#find_package(GSL REQUIRED)
#use the ROOT things
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#list(APPEND CMAKE_PREFIX_PATH $ENV{CLHEP_ROOT_DIR})
find_package(ROOT REQUIRED COMPONENTS RIO Tree MathMore)
include(${ROOT_USE_FILE})
#add_definitions(${ROOT_CXX_FLAGS})

if(DEFINED ENV{T576_SOURCE_DIR})
  set(T576_INCLUDE_DIR $ENV{T576_SOURCE_DIR}/include)
else()
  set(T576_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
endif()

if(DEFINED ENV{T576_INSTALL_DIR})
  set(T576_INSTALL_DIR $ENV{T576_INSTALL_DIR})
  set(LIB_INSTALL_DIR $ENV{T576_INSTALL_DIR}/lib)
else()
  set(T576_INSTALL_DIR /usr/local/)
  set(LIB_INSTALL_DIR /usr/local/lib)
endif()


include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include ${CLHEP_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${T576_INCLUDE_DIR} ${T576_INCLUDE_DIR}/../)

#----------------------------------------------------------------------------


file(GLOB SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cc ${CLHEP_INCLUDE_DIR}/Vector/*.icc ${CLHEP_INCLUDE_DIR}/Units/*.icc)
file(GLOB HEADER_FILES ${T576_INCLUDE_DIR}/*.hh ${CLHEP_INCLUDE_DIR}/Vector/*.h ${CLHEP_INCLUDE_DIR}/Units/*.h)# ${GSL_INCLUDE_DIRS}/gsl/*.h)


#make the ROOT dictionaries and shared object library for reading the event tree

ROOT_GENERATE_DICTIONARY(G__t576 ${HEADER_FILES} LINKDEF ${T576_INCLUDE_DIR}/../LinkDef.h)


add_library(t576 SHARED ${SOURCE_FILES} G__t576.cxx)

target_link_libraries(t576 ${ROOT_LIBRARIES} CLHEP::Vector cnpy)


#----------------------------------------------------------------------------


install (FILES ${T576_INCLUDE_DIR}/T576Event.hh ${T576_INCLUDE_DIR}/TUtil.hh ${T576_INCLUDE_DIR}/TUtilGraph.hh DESTINATION ${T576_INSTALL_DIR}/include/t576)
install (FILES ${PROJECT_SOURCE_DIR}/src/T576Event.cc DESTINATION ${T576_INSTALL_DIR}/src/t576)

install (FILES ${PROJECT_SOURCE_DIR}/share/runLog.root DESTINATION ${T576_INSTALL_DIR}/share/t576)

install (FILES ${PROJECT_BINARY_DIR}/libt576_rdict.pcm DESTINATION ${LIB_INSTALL_DIR})
install (TARGETS t576 LIBRARY DESTINATION ${LIB_INSTALL_DIR})
#install(TARGETS loadEvents DESTINATION ${T576_INSTALL_DIR}/share/t576)







