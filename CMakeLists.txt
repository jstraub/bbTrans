# Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu>
# Licensed under the MIT license. See the license file LICENSE.

project(bbTrans)
cmake_minimum_required(VERSION 2.8.1)
set(POD_NAME bbTrans)
include(cmake/pods.cmake)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_COLOR_MAKEFILE ON)
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug) # Release; debug; ...
  set(CMAKE_BUILD_TYPE Release) # Release; debug; ...
ENDIF()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

find_package (Eigen3 REQUIRED)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_RELEASE "-g -Wall -Wno-deprecated -std=c++11 -O2 -finline-functions -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "-g -Wall -Wno-deprecated -O2 -finline-functions")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wno-deprecated -std=c++11")

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  )

add_library(bbTrans SHARED
  ./src/tetrahedron.cpp
  ./src/node.cpp
  ./src/node_S3.cpp
  ./src/node_Lin.cpp
  ./src/node_TpS3.cpp
  ./src/node_AA.cpp
  ./src/node_R3.cpp
  ./src/box.cpp
  ./src/combinations.cpp
  ./src/s3_tessellation.cpp
  ./src/numeric_helpers.cpp
  ./src/lower_bound_S3.cpp
  ./src/upper_bound_indep_S3.cpp
  ./src/upper_bound_convex_S3.cpp
  ./src/lower_bound_R3.cpp
  ./src/upper_bound_indep_R3.cpp
  ./src/upper_bound_convex_R3.cpp
  )

add_executable(testCombinations ./test/combinations.cpp)
target_link_libraries(testCombinations bbTrans)

add_executable(testS3Tessellation ./test/s3_tessellation.cpp)
target_link_libraries(testS3Tessellation bbTrans)

add_executable(testBoundsS3 ./test/bounds_S3.cpp)
target_link_libraries(testBoundsS3 bbTrans)

add_executable(testBoundsR3 ./test/bounds_R3.cpp)
target_link_libraries(testBoundsR3 bbTrans)

add_executable(testBranchAndBoundS3 ./test/branch_and_bound_S3.cpp)
target_link_libraries(testBranchAndBoundS3 bbTrans)

add_executable(testBranchAndBoundTpS3 ./test/branch_and_bound_TpS3.cpp)
target_link_libraries(testBranchAndBoundTpS3 bbTrans)

add_executable(testBranchAndBoundAA ./test/branch_and_bound_AA.cpp)
target_link_libraries(testBranchAndBoundAA bbTrans)

add_executable(testBranchAndBoundR3 ./test/branch_and_bound_R3.cpp)
target_link_libraries(testBranchAndBoundR3 bbTrans)

add_executable(testBranchAndBoundT3 ./test/branch_and_bound_T3.cpp)
target_link_libraries(testBranchAndBoundT3 bbTrans)

add_executable(testMinMaxInBox ./test/min_max_in_box.cpp)
target_link_libraries(testMinMaxInBox bbTrans)

add_executable(testMinMaxInTetrahedron ./test/min_max_in_tetrahedron.cpp)
target_link_libraries(testMinMaxInTetrahedron bbTrans)

pods_install_libraries(bbTrans)
pods_install_headers(
  ./include/bbTrans/s3_tessellation.h
  ./include/bbTrans/tetrahedron.h
  ./include/bbTrans/box.h
  ./include/bbTrans/combinations.h
  ./include/bbTrans/bound.h
  ./include/bbTrans/node.h
  ./include/bbTrans/node_impl.h
  ./include/bbTrans/node_S3.h
  ./include/bbTrans/node_R3.h
  ./include/bbTrans/node_Lin.h
  ./include/bbTrans/node_TpS3.h
  ./include/bbTrans/node_AA.h
  ./include/bbTrans/branch_and_bound.h
  ./include/bbTrans/branch_and_bound_impl.h
  ./include/bbTrans/numeric_helpers.h
  ./include/bbTrans/vmf.h
  ./include/bbTrans/vmf_impl.h
  ./include/bbTrans/vmf_mm.h
  ./include/bbTrans/vmf_mm_impl.h
  ./include/bbTrans/normal.h
  ./include/bbTrans/normal_impl.h
  ./include/bbTrans/lower_bound_R3.h
  ./include/bbTrans/upper_bound_indep_R3.h
  ./include/bbTrans/upper_bound_convex_R3.h
  ./include/bbTrans/lower_bound_S3.h
  ./include/bbTrans/upper_bound_indep_S3.h
  ./include/bbTrans/upper_bound_convex_S3.h
  ./include/bbTrans/lower_bound_Lin.h
  ./include/bbTrans/lower_bound_Lin_impl.hpp
  ./include/bbTrans/upper_bound_Lin.h
  ./include/bbTrans/upper_bound_Lin_impl.hpp
  DESTINATION bbTrans)
