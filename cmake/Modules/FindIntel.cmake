# - Find the Intel MKL/IPP libraries, and Intel OpenMP library if necessary
# Modified from Armadillo's ARMA_FindMKL.cmake
# This module defines
#  MKL_INCLUDE_DIR, the directory for the MKL headers
#  MKL_LIB_DIR, the directory for the MKL library files
#  IPP_INCLUDE_DIR, the directory for the IPP headers
#  IPP_LIB_DIR, the directory for the IPP library files
#  INTEL_COMPILER_LIB_DIR, the directory for the MKL compiler library files
#  INTEL_LIBRARIES, the libraries needed to use Intel's implementation of BLAS & LAPACK.
#  INTEL_FOUND, If false, do not try to use.

# NEW ----
# 
# - If Intel compiler/MKL/IPP are installed as part of Intel oneAPI toolkit,
#	directory structure is different to older installs (at least on linux)
#	
#	2 extra directory layers added:
#   - /oneapi inside /intel
#	  so MKL etc. placed in /opt/intel/oneapi instead of /opt/intel
#	- subdirectories for version (e.g. 2022.0.2) placed in /oneapi/mkl, /oneapi/ipp etc.
#	  (symlink /oneapi/latest automatically points to these)
#   
# 	INTEL_ROOT_PATH now includes /oneapi if present (Mac & Linux)
#	(and will ignore any other /opt/intel subdirectories)
# 	
#	MKL_INCLUDE_DIR, IPP_INCLUDE_DIR now include /latest if present
# 	
#	MKL_LIB_DIR, IPP_LIB_DIR now include /latest if present
#	(not on Apple - unsure if this needs changing)
#
# 	(Linux only) INTEL_COMPILER_LIB_DIR now includes /latest if present, and checks additional path
#	
#	Above changes for 64-bit only, unsure if should be changed for 32-bit
#
# - If not installed by root, Intel compiler, MKL & IPP install in home directory
# 	INTEL_ROOT_PATH now looks in $ENV{HOME} after looking in /opt
# xxx ----


if(WIN32)

  if(NOT DEFINED ENV{INTEL_ROOT_PATH})
    set(INTEL_ROOT_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows" CACHE PATH "Where the MKL are stored")
  endif()

else()

  if(DEFINED ENV{INTEL_ROOT_PATH})
    set(INTEL_ROOT_PATH "$ENV{INTEL_ROOT_PATH}" CACHE PATH "Where the MKL are stored")
  else()
    if(EXISTS /opt/intel/oneapi)
    	set(INTEL_ROOT_PATH "/opt/intel/oneapi" CACHE PATH "Where the MKL are stored")
    elseif(EXISTS /opt/intel)
    	set(INTEL_ROOT_PATH "/opt/intel" CACHE PATH "Where the MKL are stored")
    elseif(EXISTS $ENV{HOME}/intel/oneapi)
    	set(INTEL_ROOT_PATH "$ENV{HOME}/intel/oneapi" CACHE PATH "Where the MKL are stored")
    else()
    	set(INTEL_ROOT_PATH "$ENV{HOME}/intel" CACHE PATH "Where the MKL are stored")
    endif()
  endif()

endif()


if(EXISTS ${INTEL_ROOT_PATH}/mkl)

  set(INTEL_FOUND TRUE)
  
  if(EXISTS ${INTEL_ROOT_PATH}/mkl/latest)
  	message("MKL is found at ${INTEL_ROOT_PATH}/mkl/latest")
  else()
  	message("MKL is found at ${INTEL_ROOT_PATH}/mkl")
  endif()

  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set( USE_MKL_64BIT On )
  else(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set( USE_MKL_64BIT Off )
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

else()

  set(INTEL_FOUND FALSE)
  message("MKL is NOT found ... ")

endif()


if(EXISTS ${INTEL_ROOT_PATH}/ipp)

  if(EXISTS ${INTEL_ROOT_PATH}/ipp/latest)
  	message("IPP is found at ${INTEL_ROOT_PATH}/ipp/latest")
  else()
  	message("IPP is found at ${INTEL_ROOT_PATH}/ipp")
  endif()

else()

  set(INTEL_FOUND FALSE)
  message("IPP is NOT found ... ")

endif()


if(INTEL_FOUND)

  if(EXISTS ${INTEL_ROOT_PATH}/mkl/latest)
    set(MKL_INCLUDE_DIR "${INTEL_ROOT_PATH}/mkl/latest/include")
  else()
  	set(MKL_INCLUDE_DIR "${INTEL_ROOT_PATH}/mkl/include")
  endif()
  
  if(EXISTS ${INTEL_ROOT_PATH}/ipp/latest)
    set(IPP_INCLUDE_DIR "${INTEL_ROOT_PATH}/ipp/latest/include")
  else()
  	set(IPP_INCLUDE_DIR "${INTEL_ROOT_PATH}/ipp/include")
  endif()
  
  add_definitions(-DUSE_MKL)

  if(USE_MKL_64BIT)

    if(APPLE)
      set(MKL_LIB_DIR "${INTEL_ROOT_PATH}/mkl/lib")
      set(IPP_LIB_DIR "${INTEL_ROOT_PATH}/ipp/lib")
      set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/lib")
    else()
    
      if(EXISTS ${INTEL_ROOT_PATH}/mkl/latest)
    	set(MKL_LIB_DIR "${INTEL_ROOT_PATH}/mkl/latest/lib/intel64")
      else()
    	set(MKL_LIB_DIR "${INTEL_ROOT_PATH}/mkl/lib/intel64")	
      endif()
        
      if(EXISTS ${INTEL_ROOT_PATH}/ipp/latest)
    	set(IPP_LIB_DIR "${INTEL_ROOT_PATH}/ipp/latest/lib/intel64")
      else()
    	set(IPP_LIB_DIR "${INTEL_ROOT_PATH}/ipp/lib/intel64")	
      endif()
            
      if(WIN32)
        set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/compiler/lib/intel64")
      else()
      	if(EXISTS ${INTEL_ROOT_PATH}/compiler/latest/linux/compiler/lib/intel64)
      		set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/compiler/latest/linux/compiler/lib/intel64")
      	else()
      		set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/lib/intel64")
      	endif()
      	
      endif()
      
    endif()

    if(USE_MKL_64BIT_LIB)
      if(APPLE)
        set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_ilp64.dylib)
      else()
        if(WIN32)
          set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_ilp64.lib)
        else()
          set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_ilp64.so)
        endif()
      endif()
    else()
      if(APPLE)
        set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_lp64.dylib)
      else(APPLE)
        if(WIN32)
          set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_lp64.lib)
        else()
          set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_lp64.so)
        endif()
      endif()
    endif()

  else()

    set(MKL_LIB_DIR "${INTEL_ROOT_PATH}/mkl/lib/ia32")
    set(IPP_LIB_DIR "${INTEL_ROOT_PATH}/ipp/lib/ia32")
    if(WIN32)
      set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/compiler/lib/ia32")
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_c.lib)
    else()
      set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/lib/ia32")
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel.so)
    endif()

  endif()

  if(APPLE)
    set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_thread.dylib)
    set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_core.dylib)
    set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${IPP_LIB_DIR}/libipps.dylib)
    set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${IPP_LIB_DIR}/libippcore.dylib)
    set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${INTEL_COMPILER_LIB_DIR}/libiomp5.dylib)
  else()
    if(WIN32)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_thread.lib)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/mkl_core.lib)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${IPP_LIB_DIR}/ipps.lib)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${IPP_LIB_DIR}/ippcore.lib)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${INTEL_COMPILER_LIB_DIR}/libiomp5md.lib)
   else()
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_thread.so)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_core.so)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${IPP_LIB_DIR}/libipps.a)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${IPP_LIB_DIR}/libippcore.a)
      set(INTEL_LIBRARIES ${INTEL_LIBRARIES} ${INTEL_COMPILER_LIB_DIR}/libiomp5.so)
    endif()
  endif()

  if(NOT INTEL_FIND_QUIETLY)
    message(STATUS "Found Intel libraries: ${INTEL_LIBRARIES}")
    message(STATUS "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}")
    message(STATUS "MKL_LIB_DIR: ${MKL_LIB_DIR}")
    message(STATUS "IPP_INCLUDE_DIR: ${IPP_INCLUDE_DIR}")
    message(STATUS "IPP_LIB_DIR: ${IPP_LIB_DIR}")
    message(STATUS "INTEL_COMPILER_LIB_DIR: ${INTEL_COMPILER_LIB_DIR}")
  endif()

else()

  if(INTEL_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find Intel MKL and/or IPP libraries")
  endif()

endif()
