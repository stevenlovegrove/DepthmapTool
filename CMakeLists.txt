##############################################################################
# CUDA configuration
PROJECT("ConvertDepthMapIntoImage")


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

# Overide with cmake -DCMAKE_BUILD_TYPE=Debug {dir}
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Release )
ENDIF()

#######################################################
## Find required dependencies

FIND_PACKAGE(CVD REQUIRED)

#######################################################

INCLUDE_DIRECTORIES(
  ${CVD_INCLUDE_DIR}
)

LINK_LIBRARIES(
  ${CVD_LIBRARY}
)

ADD_EXECUTABLE(
  ConvertDepthMapIntoImage
  createdepthmapimage.cpp
)
