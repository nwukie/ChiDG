
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
find_library (TECPLOT_LIBRARIES
        NAMES libtecio.a libtecio
        HINTS ${TEC_ROOT}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)


# SHOULD PROBABLY REDO THIS
#function (set_tecplot_install_dir)
#    set (TEC_HOME $ENV{TEC_ROOT})
#    if (TEC_HOME)
#        set (TEC_HOME $ENV{TEC_ROOT})
#    else()
#        message (FATAL_ERROR "TEC_ROOT does not exist. Please set it to the location of the Tecplot installation directory.")
#    endif ()
#
#endfunction ()
#
#
#set_tecplot_install_dir ()
#
#
#find_path (TECPLOT_INCLUDE_DIR
#           "tecio.f90"
#           HINTS $ENV{TEC_ROOT}
#           PATH_SUFFIXES include
#           NO_DEFAULT_PATH)
#find_library (TECPLOT_LIBRARIES
#              NAMES libtecio.a libtecio
#              HINTS $ENV{TEC_ROOT}
#              PATH_SUFFIXES lib
#              NO_DEFAULT_PATH)
#
#message(STATUS "TEC_ROOT: " $ENV{TEC_ROOT})
#message(STATUS "TEC_HOME: " ${TEC_HOME})
#
#set(TECPLOT_INCLUDE_DIR $ENV{TEC_ROOT}/include)
#set(TECPLOT_LIBRARIES $ENV{TEC_ROOT}/lib/libtecio.a)




