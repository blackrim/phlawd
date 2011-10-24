cmake_minimum_required(VERSION 2.8)

project(phlawd CXX)

include(FindPkgConfig)

# Borrowed from the BPP code.
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   libraries: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

IMPROVED_FIND_LIBRARY(LIBS bpp-seq Seq/Sequence.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-numcalc NumCalc/NumTools.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-utils Utils/Clonable.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Phyl/Node.h)
pkg_search_module(SQLITE3 REQUIRED sqlite3)

add_subdirectory(src)
