# - Find NetCDF
# Find the native NetCDF includes and library
#
#  netCDF_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  netCDF_LIBRARY      - user modifiable choice of where netcdf libraries are
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  netCDF_CXX         - require the C++ interface and link the C++ library
#  netCDF_F77         - require the F77 interface and link the fortran library
#  netCDF_F90         - require the F90 interface and link the fortran library
#
# Or equivalently by calling FindNetCDF with a COMPONENTS argument containing one or
# more of "CXX;F77;F90".
#
# When interfaces are requested the user has access to interface specific hints:
#
#  netCDF_${LANG}_INCLUDE_DIR - where to search for interface header files
#  netCDF_${LANG}_LIBRARY     - where to search for interface libraries
#
# This module returns these variables for the rest of the project to use.
#
#  netCDF_FOUND          - True if NetCDF found including required interfaces (see below)
#  netCDF_LIBRARIES      - All netcdf related libraries.
#  netCDF_INCLUDE_DIRS   - All directories to include.
#  netCDF_HAS_INTERFACES - Whether requested interfaces were found or not.
#  netCDF_${LANG}_INCLUDE_DIRS/netCDF_${LANG}_LIBRARIES - C/C++/F70/F90 only interface
#
# Normal usage would be:
#  set (netCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_everthing ${netCDF_LIBRARIES})
#  target_link_libraries (only_uses_f90 ${netCDF_F90_LIBRARIES})

#search starting from user editable cache var
if (netCDF_INCLUDE_DIR AND netCDF_LIBRARY)
  # Already in cache, be silent
  set (netCDF_FIND_QUIETLY TRUE)
endif ()

set(USE_DEFAULT_PATHS "NO_DEFAULT_PATH")
if(netCDF_USE_DEFAULT_PATHS)
  set(USE_DEFAULT_PATHS "")
endif()

find_path (netCDF_INCLUDE_DIR netcdf.h
  HINTS "${netCDF_DIR}/../../../include")
mark_as_advanced (netCDF_INCLUDE_DIR)
set (netCDF_C_INCLUDE_DIRS ${netCDF_INCLUDE_DIR})

find_library (netCDF_LIBRARY NAMES netcdf
  HINTS "${netCDF_DIR}/../..")
mark_as_advanced (netCDF_LIBRARY)

set (netCDF_C_LIBRARIES ${netCDF_LIBRARY})

#start finding requested language components
set (NetCDF_libs "")
set (NetCDF_includes "${netCDF_INCLUDE_DIR}")

get_filename_component (NetCDF_lib_dirs "${netCDF_LIBRARY}" PATH)
set (netCDF_HAS_INTERFACES "YES") # will be set to NO if we're missing any interfaces

macro (NetCDF_check_interface lang header libs)
  if (netCDF_${lang})
    #search starting from user modifiable cache var
    find_path (netCDF_${lang}_INCLUDE_DIR NAMES ${header}
      HINTS "${netCDF_INCLUDE_DIR}"
      HINTS "${netCDF_${lang}_ROOT}/include"
      ${USE_DEFAULT_PATHS})

    find_library (netCDF_${lang}_LIBRARY NAMES ${libs}
      HINTS "${NetCDF_lib_dirs}"
      HINTS "${netCDF_${lang}_ROOT}/lib"
      ${USE_DEFAULT_PATHS})

    mark_as_advanced (netCDF_${lang}_INCLUDE_DIR netCDF_${lang}_LIBRARY)

    #export to internal varS that rest of project can use directly
    set (netCDF_${lang}_LIBRARIES ${netCDF_${lang}_LIBRARY})
    set (netCDF_${lang}_INCLUDE_DIRS ${netCDF_${lang}_INCLUDE_DIR})

    if (netCDF_${lang}_INCLUDE_DIR AND netCDF_${lang}_LIBRARY)
      list (APPEND NetCDF_libs ${netCDF_${lang}_LIBRARY})
      list (APPEND NetCDF_includes ${netCDF_${lang}_INCLUDE_DIR})
    else ()
      set (netCDF_HAS_INTERFACES "NO")
      message (STATUS "Failed to find NetCDF interface for ${lang}")
    endif ()
  endif ()
endmacro ()

list (FIND NetCDF_FIND_COMPONENTS "CXX" _nextcomp)
if (_nextcomp GREATER -1)
  set (netCDF_CXX 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F77" _nextcomp)
if (_nextcomp GREATER -1)
  set (netCDF_F77 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F90" _nextcomp)
if (_nextcomp GREATER -1)
  set (netCDF_F90 1)
endif ()
NetCDF_check_interface (CXX netcdfcpp.h netcdf_c++)
NetCDF_check_interface (F77 netcdf.inc  netcdff)
NetCDF_check_interface (F90 netcdf.mod  netcdff)

#export accumulated results to internal varS that rest of project can depend on
list (APPEND NetCDF_libs "${netCDF_C_LIBRARIES}")
set (netCDF_LIBRARIES ${NetCDF_libs})
set (netCDF_INCLUDE_DIRS ${NetCDF_includes})

# handle the QUIETLY and REQUIRED arguments and set netCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF
  DEFAULT_MSG netCDF_LIBRARIES netCDF_INCLUDE_DIRS netCDF_HAS_INTERFACES)
