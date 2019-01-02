# $Id: CMakeLists.txt 100923 2016-11-03 10:50:34Z gcosmo $

#



#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(t576)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
unset(CMAKE_MODULE_PATH)



find_package(ZLIB REQUIRED)



#find_package(GSL REQUIRED)
#use the ROOT things
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

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




include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${T576_INCLUDE_DIR} ${T576_INCLUDE_DIR}/../)

#----------------------------------------------------------------------------




file(GLOB SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cpp )
file(GLOB HEADER_FILES ${T576_INCLUDE_DIR}/*.h ${T576_INCLUDE_DIR}/*.hh)# ${GSL_INCLUDE_DIRS}/gsl/*.h)


#make the ROOT dictionaries and shared object library for reading the event tree

ROOT_GENERATE_DICTIONARY(G__t576 ${HEADER_FILES} LINKDEF ${T576_INCLUDE_DIR}/../LinkDef.h)


add_library(t576 SHARED ${SOURCE_FILES} G__t576.cxx)


target_link_libraries(t576 ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})# cnpy)
add_executable(svdExample ${PROJECT_SOURCE_DIR}/example/svdExample.cc)
target_link_libraries(svdExample ${ROOT_LIBRARIES} t576)
#----------------------------------------------------------------------------


install (FILES ${T576_INCLUDE_DIR}/T576Event.hh ${T576_INCLUDE_DIR}/TUtil.hh ${T576_INCLUDE_DIR}/cnpy.h ${T576_INCLUDE_DIR}/TUtilGraph.hh DESTINATION ${T576_INSTALL_DIR}/include/t576)
install (FILES ${PROJECT_SOURCE_DIR}/src/T576Event.cc DESTINATION ${T576_INSTALL_DIR}/src/t576)

install (FILES ${PROJECT_SOURCE_DIR}/share/runLog.root DESTINATION ${T576_INSTALL_DIR}/share/t576)

install (FILES ${PROJECT_BINARY_DIR}/libt576_rdict.pcm DESTINATION ${LIB_INSTALL_DIR})
install (TARGETS t576 LIBRARY DESTINATION ${LIB_INSTALL_DIR})
install(TARGETS svdExample DESTINATION ${T576_INSTALL_DIR}/share/t576/example)







