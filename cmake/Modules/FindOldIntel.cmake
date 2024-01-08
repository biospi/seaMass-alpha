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


if(WIN32)

  if(NOT DEFINED ENV{INTEL_ROOT_PATH})
    set(INTEL_ROOT_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows" CACHE PATH "Where the MKL are stored")
  endif()

else()

  if(DEFINED ENV{INTEL_ROOT_PATH})
    set(INTEL_ROOT_PATH "$ENV{INTEL_ROOT_PATH}" CACHE PATH "Where the MKL are stored")
  else()
    set(INTEL_ROOT_PATH "/opt/intel" CACHE PATH "Where the MKL are stored")
  endif()

endif()


if(EXISTS ${INTEL_ROOT_PATH}/mkl)

  set(INTEL_FOUND TRUE)
  message("MKL is found at ${INTEL_ROOT_PATH}/mkl")

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

  message("IPP is found at ${INTEL_ROOT_PATH}/ipp")

else()

  set(INTEL_FOUND FALSE)
  message("IPP is NOT found ... ")

endif()


if(INTEL_FOUND)

  set(MKL_INCLUDE_DIR "${INTEL_ROOT_PATH}/mkl/include")
  set(IPP_INCLUDE_DIR "${INTEL_ROOT_PATH}/ipp/include")
  add_definitions(-DUSE_MKL)

  if(USE_MKL_64BIT)

    if(APPLE)
      set(MKL_LIB_DIR "${INTEL_ROOT_PATH}/mkl/lib")
      set(IPP_LIB_DIR "${INTEL_ROOT_PATH}/ipp/lib")
      set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/lib")
    else()
      set(MKL_LIB_DIR "${INTEL_ROOT_PATH}/mkl/lib/intel64")
      set(IPP_LIB_DIR "${INTEL_ROOT_PATH}/ipp/lib/intel64")
      if(WIN32)
        set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/compiler/lib/intel64")
      else()
        set(INTEL_COMPILER_LIB_DIR "${INTEL_ROOT_PATH}/lib/intel64")
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
