# - Find the MKL libraries, and Intel OpenMP library if necessary
# Modified from Armadillo's ARMA_FindMKL.cmake
# This module defines
#  MKL_INCLUDE_DIR, the directory for the MKL headers
#  MKL_LIB_DIR, the directory for the MKL library files
#  MKL_COMPILER_LIB_DIR, the directory for the MKL compiler library files
#  MKL_LIBRARIES, the libraries needed to use Intel's implementation of BLAS & LAPACK.
#  MKL_FOUND, If false, do not try to use MKL; if true, the macro definition USE_MKL is added.

# Set the include path
# TODO: what if MKL is not installed in /opt/intel/mkl?
# try to find at /opt/intel/mkl
# in windows, try to find MKL at C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl


if(WIN32)

  if(NOT DEFINED ENV{MKLROOT_PATH})
    set(MKLROOT_PATH "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows" CACHE PATH "Where the MKL are stored")
  endif(NOT DEFINED ENV{MKLROOT_PATH}) 

else(WIN32)

  if(DEFINED ENV{MKLROOT_PATH})
    set(MKLROOT_PATH "$ENV{MKLROOT_PATH}" CACHE PATH "Where the MKL are stored")
  else(DEFINED ENV{MKLROOT_PATH})
    set(MKLROOT_PATH "/opt/intel" CACHE PATH "Where the MKL are stored")
  endif(DEFINED ENV{MKLROOT_PATH})

endif(WIN32)


if(EXISTS ${MKLROOT_PATH}/mkl)

  set(MKL_FOUND TRUE)
  message("MKL is found at ${MKLROOT_PATH}/mkl")

  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set( USE_MKL_64BIT On )
  else(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set( USE_MKL_64BIT Off )
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

else(EXISTS ${MKLROOT_PATH}/mkl)

  set(MKL_FOUND FALSE)
  message("MKL is NOT found ... ")

endif(EXISTS ${MKLROOT_PATH}/mkl)


if(MKL_FOUND)

  set(MKL_INCLUDE_DIR "${MKLROOT_PATH}/mkl/include")
  add_definitions(-DUSE_MKL)

  if(USE_MKL_64BIT)

    if(APPLE)
      set(MKL_LIB_DIR "${MKLROOT_PATH}/mkl/lib")
      set(MKL_COMPILER_LIB_DIR "${MKLROOT_PATH}/compiler/lib")
      set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT_PATH}/lib")
    else(APPLE)
      set(MKL_LIB_DIR "${MKLROOT_PATH}/mkl/lib/intel64")
      set(MKL_COMPILER_LIB_DIR "${MKLROOT_PATH}/compiler/lib/intel64")
      set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT_PATH}/lib/intel64")
    endif(APPLE)

    if(USE_MKL_64BIT_LIB)
      if(APPLE)
        set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_ilp64.dylib)
      else(APPLE)
        if(WIN32)
          set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_ilp64.lib)
        else(WIN32)
          set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_ilp64.so)
        endif(WIN32)
      endif(APPLE)
    else(USE_MKL_64BIT_LIB)
      if(APPLE)
        set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_lp64.dylib)
      else(APPLE)
        if(WIN32)
          set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_lp64.lib)
        else(WIN32)
          set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_lp64.so)
        endif(WIN32)
      endif(APPLE)
    endif(USE_MKL_64BIT_LIB)

  else(USE_MKL_64BIT)

    set(MKL_LIB_DIR "${MKLROOT_PATH}/mkl/lib/ia32")
    set(MKL_COMPILER_LIB_DIR "${MKLROOT_PATH}/compiler/lib/ia32")
    set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT_PATH}/lib/ia32")
    if (WIN32)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_c.lib)
    else (WIN32)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel.so)
    endif (WIN32)

  endif(USE_MKL_64BIT)

  if(APPLE)
    set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_intel_thread.dylib)
    set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_core.dylib)
    if(OPENMP_FOUND)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKLROOT_PATH}/lib/libiomp5.dylib)
    endif(OPENMP_FOUND)
  else(APPLE)
    if(WIN32)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/mkl_intel_thread.lib)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/mkl_core.lib)
    else (WIN32)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_gnu_thread.so)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_LIB_DIR}/libmkl_core.so)
    endif(WIN32) 
  endif(APPLE)

  if(NOT MKL_FIND_QUIETLY)
    message(STATUS "Found MKL libraries: ${MKL_LIBRARIES}")
    message(STATUS "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}")
    message(STATUS "MKL_LIB_DIR: ${MKL_LIB_DIR}")
    message(STATUS "MKL_COMPILER_LIB_DIR: ${MKL_COMPILER_LIB_DIR}")
  endif(NOT MKL_FIND_QUIETLY)

else(MKL_FOUND)

  if(MKL_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find MKL libraries")
  endif(MKL_FIND_REQUIRED)

endif(MKL_FOUND)
