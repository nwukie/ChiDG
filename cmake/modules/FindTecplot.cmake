
# Find Tecplot
# Finds the Tecplot libraries and include directory needed for an add-on.
# The critical variable is TEC_ROOT. Once that is set, this
# module will find with include directory and libraries within the
# install directory.
#
# Typical usage:
#   find_package (Tecplot REQUIRED)
#   include_directories ("${TECPLOT_INCLUDE_DIR}")


find_path (TECPLOT_INCLUDE_DIR
           "tecio.f90"
           HINTS ${TEC_ROOT}
           PATH_SUFFIXES include
           NO_DEFAULT_PATH)
find_library (TECIO_LIBRARY
              NAMES libtecio.a libtecio.so libtecio.dylib libtecio
              HINTS ${TEC_ROOT}
              PATH_SUFFIXES lib
              NO_DEFAULT_PATH)


#
# Append the TECPLOT_LIBRARIES with the Threads and C++ libraries
#
set(TECPLOT_LIBRARIES ${TECIO_LIBRARY} "stdc++" "pthread")


