# Find the pugixml XML parsing library.
#
# Sets the usual variables expected for find_package scripts:
#
# PUGIXML_INCLUDE_DIR - header location
# PUGIXML_LIBRARIES - library to link against
# PUGIXML_FOUND - true if pugixml was found.
#
# If PugiXML is located in another location the following Environment variable:
# export PUGIXML_ROOT_PATH=/location/of/pugiXML

if(DEFINED ENV{PUGIXML_ROOT_PATH})
    set(PUGIXML_ROOT_PATH "$ENV{PUGIXML_ROOT_PATH}" CACHE PATH "Where the PugiXML is stored")
    message("PugiXML defined at ${PUGIXML_ROOT_PATH}")
    find_path (PugiXML_INCLUDE_DIR NAMES pugixml.hpp PATHS
        ${PUGIXML_ROOT_PATH}/include
    )
    find_library (PugiXML_LIBRARY NAMES pugixml PATHS
        ${PUGIXML_ROOT_PATH}/lib
    )
else()
    message("PugiXML Normal Search")
    file(GLOB PugiXML_INCLUDE_SEARCH_PATH "/usr/local/include/pugixml-*")
    find_path (PugiXML_INCLUDE_DIR NAMES pugixml.hpp PATHS
        ${PugiXML_DIR}/../../../include
        ${PugiXML_INCLUDE_SEARCH_PATH}
    )
    file(GLOB PugiXML_LIB_SEARCH_PATH "/usr/local/lib/pugixml-*")
    find_library (PugiXML_LIBRARY NAMES pugixml PATHS
        ${PugiXML_DIR}/../..
        /usr/local/lib/pugixml-1.8
        ${PugiXML_LIB_SEARCH_PATH}
    )
endif()


# Support the REQUIRED and QUIET arguments, and set PUGIXML_FOUND if found.
include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (PugiXML FOUND_VAR PugiXML_FOUND REQUIRED_VARS PugiXML_INCLUDE_DIR PugiXML_LIBRARY)

if (PugiXML_FOUND)
    set (PugiXML_LIBRARIES ${PugiXML_LIBRARY})
    if (NOT PugiXML_FIND_QUIETLY)
        message (STATUS "PugiXML include = ${PugiXML_INCLUDE_DIR}")
        message (STATUS "PugiXML library = ${PugiXML_LIBRARY}")
    endif ()
else (PugiXML_FOUND)
    message (STATUS "No PugiXML found")
endif(PugiXML_FOUND)

mark_as_advanced (PugiXML_LIBRARY PugiXML_INCLUDE_DIR)
