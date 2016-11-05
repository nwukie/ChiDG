# Find METIS
#
# Finds the METIS libraries and include directory.
#
#
# Typical usage:
#   find_package (METIS REQUIRED)
#   include_directories (${METIS_INCLUDE_DIR})
#
#


# Set path hints from environment variable, METIS_ROOT if it is set
set(METIS_HINT_ROOT $ENV{METIS_ROOT})


# DEPRECATED APPROACH: Set path hints from environment variables if they are set
set(METIS_HINT_LIBRARY_DIR $ENV{METIS_LIBRARY_DIR})
set(METIS_HINT_INCLUDE_DIR $ENV{METIS_INCLUDE_DIR})


# Find include file and library
find_path(METIS_INCLUDE_DIR 
          metis.h 
          HINTS ${METIS_HINT_ROOT} ${METIS_HINT_INCLUDE_DIR})

find_library(METIS_LIBRARIES 
             metis 
             HINTS ${METIS_HINT_ROOT} ${METIS_HINT_LIBRARY_DIR})

# Set METIS_FOUND if conditions are met
if(METIS_INCLUDE_DIR)
    if(METIS_LIBRARIES)
        set(METIS_FOUND TRUE) 
    endif(METIS_LIBRARIES)
endif(METIS_INCLUDE_DIR)



