cmake_minimum_required(VERSION 2.8)

set(project_policies
  CMP0001 # NEW: CMAKE_BACKWARDS_COMPATIBILITY should no longer be used.
  CMP0002 # NEW: Logical target names must be globally unique.
  CMP0003 # NEW: Libraries linked via full path no longer produce linker search paths.
  CMP0004 # NEW: Libraries linked may NOT have leading or trailing whitespace.
  CMP0005 # NEW: Preprocessor definition values are now escaped automatically.
  CMP0006 # NEW: Installing MACOSX_BUNDLE targets requires a BUNDLE DESTINATION.
  CMP0007 # NEW: List command no longer ignores empty elements.
  CMP0008 # NEW: Libraries linked by full-path must have a valid library file name.
  CMP0009 # NEW: FILE GLOB_RECURSE calls should not follow symlinks by default.
  CMP0010 # NEW: Bad variable reference syntax is an error.
  CMP0011 # NEW: Included scripts do automatic cmake_policy PUSH and POP.
  CMP0012 # NEW: if() recognizes numbers and boolean constants.
  CMP0013 # NEW: Duplicate binary directories are not allowed.
  CMP0014 # NEW: Input directories must have CMakeLists.txt
  CMP0042
  CMP0048)

foreach(policy ${project_policies})
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

cmake_policy(SET CMP0048 NEW)

PROJECT(cmic_dev_jc
  VERSION 1.0)

set(GIF_VERSION_MAJOR 3)
set(GIF_VERSION_MINOR 0)
set(GIF_VERSION_PATCH 0)
include(ExternalProject)


set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY lib)

set(LIBRARY_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(ARCHIVE_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(RUNTIME_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})

set(LIBRARY_INSTALL_OUTPUT_PATH    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(ARCHIVE_INSTALL_OUTPUT_PATH    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(RUNTIME_INSTALL_OUTPUT_PATH    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})
set(EXECUTABLE_INSTALL_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})

set(${PROJECT_NAME}_VERSION_MAJOR 1 CACHE STRING "${PROJECT_NAME} major version number.")
set(${PROJECT_NAME}_VERSION_MINOR 0 CACHE STRING "${PROJECT_NAME} minor version number.")
set(${PROJECT_NAME}_VERSION_BUILD 0 CACHE STRING "${PROJECT_NAME} build version number.")
set(${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_BUILD}")

mark_as_advanced(${PROJECT_NAME}_VERSION_MAJOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_MINOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_BUILD)

#-----------------------------------------------------------------------------
# Add compiler flags ${PROJECT_NAME} needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.
# [leukoquant patch] The vendored NiftyReg's RNifti headers use C++17
# features (std::variant, etc.) since its post-#130 refactor -- bumped from
# 11 to match NiftyReg's own CMAKE_CXX_STANDARD requirement.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(${PROJECT_NAME}_REQUIRED_C_FLAGS)
set(${PROJECT_NAME}_REQUIRED_CXX_FLAGS)
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${${PROJECT_NAME}_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${PROJECT_NAME}_REQUIRED_CXX_FLAGS}")

set(CMAKE_MODULE_PATH
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/CMake/CMakeExternals
  ${CMAKE_MODULE_PATH}
  )

set(${PROJECT_NAME}_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/Utilities
  ${PROJECT_SOURCE_DIR}/Utilities/rapidxml
  ${PROJECT_SOURCE_DIR}/Code
  ${PROJECT_SOURCE_DIR}/Code/Libraries
  ${PROJECT_SOURCE_DIR}/Code/Applications
)

include_directories(
  ${${PROJECT_NAME}_INCLUDE_DIRS}
)


option(BUILD_SHARED_LIBS "Build with shared libraries." OFF)
if(BUILD_SHARED_LIBS)
  set(${PROJECT_NAME}_BUILD_TYPE
    SHARED
    )
endif()

find_package(NIFTYREG REQUIRED)
find_package(NIFTYSEG REQUIRED)

#-----------------------------------------------------------------------------
# Search for Eigen. Download and install if not found
#-----------------------------------------------------------------------------
option(USE_SYSTEM_EIGEN "Use an already installed version of the Eigen library" OFF)
if(USE_SYSTEM_EIGEN)
  find_package(EIGEN REQUIRED)
else()
  set(${PROJECT_NAME}_VERSION_EIGEN "ffa86ffb5570" CACHE STRING "Version of EIGEN" FORCE)
  set(${PROJECT_NAME}_MD5_SUM_EIGEN a77de843845c4805b5fc18f64c74a2e6 CACHE STRING "MD5 check sum for EIGEN" FORCE)
  set(${PROJECT_NAME}_LOCATION_EIGEN 
    "https://gitlab.com/libeigen/eigen/-/archive/3.3-beta1/eigen-3.3-beta1.tar.gz"
    CACHE STRING "Location of Eigen" FORCE)
  mark_as_advanced(FORCE ${PROJECT_NAME}_VERSION_EIGEN)
  mark_as_advanced(FORCE ${PROJECT_NAME}_MD5_SUM_EIGEN)
  mark_as_advanced(FORCE ${PROJECT_NAME}_LOCATION_EIGEN)
  ExternalProject_Add(Eigen
    URL ${${PROJECT_NAME}_LOCATION_EIGEN}
    URL_MD5 ${${PROJECT_NAME}_MD5_SUM_EIGEN}
    PREFIX ${PROJECT_BINARY_DIR}/Eigen
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/Eigen/download
    SOURCE_DIR ${PROJECT_BINARY_DIR}/Eigen/source
    STAMP_DIR ${PROJECT_BINARY_DIR}/Eigen/stamps
    TMP_DIR ${PROJECT_BINARY_DIR}/Eigen/tmp
    BINARY_DIR ${PROJECT_BINARY_DIR}/Eigen/build
    CONFIGURE_COMMAND ""
    UPDATE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
  set(EIGEN_INCLUDE_DIR ${PROJECT_BINARY_DIR}/Eigen/source)
  install(DIRECTORY ${EIGEN_INCLUDE_DIR}/Eigen DESTINATION include)
endif()


include_directories(${EIGEN_INCLUDE_DIR})

if(NIFTYREG_FOUND AND NIFTYSEG_FOUND)

  include_directories(
    ${${PROJECT_NAME}_INCLUDE_DIRS}
    ${NIFTYREG_INCLUDE_DIR}
    ${NIFTYSEG_INCLUDE_DIR}
    )  

  set(CMAKE_LIBRARY_PATH 
    ${CMAKE_LIBRARY_PATH}
    ${NIFTYREG_LIBRARY_DIR}
    ${NIFTYSEG_LIBRARY_DIR}
    )

  link_directories(
    ${CMAKE_LIBRARY_PATH}
    )

  option(USE_OPENMP "To use openMP for multi-CPU processing" OFF)
  if(USE_OPENMP)
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
      set(USE_OPENMP OFF CACHE BOOL "To use openMP for multi-CPU processing" FORCE)
      message(WARNING "OpenMP does not appear to be supported by your compiler, forcing USE_OPENMP to OFF")
    else()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
  endif()

  option(USE_SSE "To use SSE for some computations" OFF)

  find_package(Git)
  if(GIT_FOUND)
          execute_process(
                  COMMAND git rev-parse HEAD
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  OUTPUT_VARIABLE GIT_COMMIT_HASH
                  RESULT_VARIABLE GIT_FAILED
                  OUTPUT_STRIP_TRAILING_WHITESPACE
          )
          if(NOT GIT_FAILED)
                   add_definitions(-D_GIT_HASH=\"${GIT_COMMIT_HASH}\")
                   message(STATUS "Found Git - the hash is added to the source.")
          else(NOT GIT_FAILED)
                   add_definitions(-D_GIT_HASH="${GIF_VERSION_MAJOR}.${GIF_VERSION_MINOR}.${GIF_VERSION_PATCH}")
          endif(NOT GIT_FAILED)
  else(GIT_FOUND)
           add_definitions(-D_GIT_HASH="${GIF_VERSION_MAJOR}.${GIF_VERSION_MINOR}.${GIF_VERSION_PATCH}")
  endif(GIT_FOUND)

  add_subdirectory(Utilities)
  add_subdirectory(Code)
  
endif()

