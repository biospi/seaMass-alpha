# - Find the Intel OpenMP library (for version 2017)
# Modified from Armadillo's ARMA_FindIntelOpenMP.cmake
# This module defines
#  IntelOpenMP_INCLUDE_DIR, the directory for the IntelOpenMP headers
#  IntelOpenMP_LIB_DIR, the directory for the IntelOpenMP library files
#  IntelOpenMP_COMPILER_LIB_DIR, the directory for the IntelOpenMP compiler library files
#  IntelOpenMP_LIBRARIES, the libraries needed to use Intel's implementation of BLAS & LAPACK.
#  IntelOpenMP_FOUND, If false, do not try to use IntelOpenMP; if true, the macro definition USE_IntelOpenMP is added.

# TODO: what if IntelOpenMP is not installed in /opt/intel/compiler?
# try to find at /opt/intel/composerxe
# in windows, try to find IntelOpenMP at C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows\compiler

if(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:vcomp")
endif()

if ( WIN32 )
  if(NOT DEFINED ENV{IntelOpenMPROOT_PATH})
    set(IntelOpenMPROOT_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler" CACHE PATH "Where the IntelOpenMP are stored")
  endif(NOT DEFINED ENV{IntelOpenMPROOT_PATH})
else ( WIN32 )
  if (DEFINED ENV{IntelOpenMPROOT_PATH})
    set(IntelOpenMPROOT_PATH "$ENV{IntelOpenMPROOT_PATH}" CACHE PATH "Where the IntelOpenMP are stored")
  else (DEFINED ENV{IntelOpenMPROOT_PATH})
    set(IntelOpenMPROOT_PATH "/opt/intel/compilers_and_libraries/linux" CACHE PATH "Where the IntelOpenMP are stored")
  endif (DEFINED ENV{IntelOpenMPROOT_PATH})
endif ( WIN32 )

if (EXISTS ${IntelOpenMPROOT_PATH})
    SET(IntelOpenMP_FOUND TRUE)
    message("IntelOpenMP is found at ${IntelOpenMPROOT_PATH}")
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( USE_IntelOpenMP_64BIT On )
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( USE_IntelOpenMP_64BIT Off )
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8)
else (EXISTS ${IntelOpenMPROOT_PATH})
    SET(IntelOpenMP_FOUND FALSE)
    message("IntelOpenMP is NOT found ... ")
endif (EXISTS ${IntelOpenMPROOT_PATH})

if (IntelOpenMP_FOUND)
    ADD_DEFINITIONS(-DUSE_IntelOpenMP)
    if ( USE_IntelOpenMP_64BIT )
        set(IntelOpenMP_LIB_DIR "${IntelOpenMPROOT_PATH}/lib/intel64")
        set(IntelOpenMP_COMPILER_LIB_DIR "${IntelOpenMPROOT_PATH}/lib/intel64")
        if (WIN32)
            set(IntelOpenMP_LIBRARIES ${IntelOpenMP_LIBRARIES} ${IntelOpenMP_LIB_DIR}/libiomp5md.lib)
        else (WIN32)
            set(IntelOpenMP_LIBRARIES ${IntelOpenMP_LIBRARIES} ${IntelOpenMP_LIB_DIR}/libiomp5.so)
        endif (WIN32)
    else ( USE_IntelOpenMP_64BIT )
        set(IntelOpenMP_LIB_DIR "${IntelOpenMPROOT_PATH}/lib/ia32")
        set(IntelOpenMP_COMPILER_LIB_DIR "${IntelOpenMPROOT_PATH}/lib/ia32")
        if ( WIN32 )
            set(IntelOpenMP_LIBRARIES ${IntelOpenMP_LIBRARIES} ${IntelOpenMP_LIB_DIR}/libiomp5md.lib)
        else ( WIN32 )
            set(IntelOpenMP_LIBRARIES ${IntelOpenMP_LIBRARIES} ${IntelOpenMP_LIB_DIR}/libiomp5.so)
        endif ( WIN32 )
    endif ( USE_IntelOpenMP_64BIT )

endif (IntelOpenMP_FOUND)

IF (IntelOpenMP_FOUND)
    IF (NOT IntelOpenMP_FIND_QUIETLY)
        MESSAGE(STATUS "Found IntelOpenMP libraries: ${IntelOpenMP_LIBRARIES}")
        MESSAGE(STATUS "IntelOpenMP_LIB_DIR: ${IntelOpenMP_LIB_DIR}")
    ENDIF (NOT IntelOpenMP_FIND_QUIETLY)
 ELSE (IntelOpenMP_FOUND)
    IF (IntelOpenMP_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find IntelOpenMP libraries")
    ENDIF (IntelOpenMP_FIND_REQUIRED)
ENDIF (IntelOpenMP_FOUND)

# MARK_AS_ADVANCED(IntelOpenMP_LIBRARY)
