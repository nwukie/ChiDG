
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


#
# Find stdc++ library used by TecIO
#
#find_library (CXX_STDLIBS
#              NAMES libstdc++ libstdc++.so libstdc++.a
#              )


#
# Find Threads library used by TecIO. 
#
#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)



#
# Append the TECPLOT_LIBRARIES with the Threads and C++ libraries
#
#set(TECPLOT_LIBRARIES ${TECPLOT_LIBRARIES} Threads::Threads "stdc++" )
set(TECPLOT_LIBRARIES ${TECPLOT_LIBRARIES} "stdc++" "pthread")


