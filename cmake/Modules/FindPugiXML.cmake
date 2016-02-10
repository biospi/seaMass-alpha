# Find the pugixml XML parsing library.
#
# Sets the usual variables expected for find_package scripts:
#
# PUGIXML_INCLUDE_DIR - header location
# PUGIXML_LIBRARIES - library to link against
# PUGIXML_FOUND - true if pugixml was found.

find_path (PugiXML_INCLUDE_DIR
           NAMES pugixml.hpp
           PATHS ${PugiXML_DIR}/../../../include)
find_library (PugiXML_LIBRARY
              NAMES pugixml
              PATHS ${PugiXML_DIR}/../..)

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
