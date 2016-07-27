project(chidg_v0.1 Fortran)
cmake_minimum_required(VERSION 3.4) 





#-----------------------------------------------------------------------------------
#                      CMake macros for adding source files
#-----------------------------------------------------------------------------------
include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)
include(cmake/macros/add_test_target.txt)
include(cmake/macros/add_test_sources.txt)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)



#-----------------------------------------------------------------------------------
#                       Enable testing with -DTESTING=True
#-----------------------------------------------------------------------------------
set(TESTING FALSE CACHE BOOL "Enable building unit/regression testing")





#-----------------------------------------------------------------------------------
#                           Get fortran compiler vendor
#-----------------------------------------------------------------------------------
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (Fortran_VENDOR "INTEL")
elseif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set (Fortran_VENDOR "GNU")
else ()
    message(SEND_ERROR "Fortran compiler not recognized")
endif()




#-----------------------------------------------------------------------------------
#                               Set compiler flags
#----------------------------------------------------------------------------------- 
if (${Fortran_VENDOR} MATCHES "INTEL")
#    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -p -g -traceback -check all -check noarg_temp_created -check bounds -fp-stack-check")
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -p -g -traceback -check all -check bounds -fp-stack-check -standard-semantics -assume nostd_mod_proc_name -heap-arrays 1600")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -fpp -assume realloc_lhs")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

    set (PROFILE_FLAGS "")
    
    #set (CMAKE_AR xiar)        # for -ipo
    #set (CMAKE_LINKER xild)    # for -ipo
    
elseif (${Fortran_VENDOR} MATCHES "GNU")
    # Note: -fcheck=all has issues with passing null pointer in tecio_interface
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp -fopenmp --coverage -ffree-form -ffree-line-length-none -fbacktrace -fstack-protector -fcheck=bounds -fmax-stack-var-size=1600 -p -g")    # Note: -fcheck=all has issues with passing null pointer in tecio_interface
    #set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -fopenmp -ffree-form -ffree-line-length-none -fexternal-blas -march=native -funroll-loops -m64 -fwhole-file")
    #set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -ffree-form -ffree-line-length-none -fexternal-blas -march=native -m64 -funroll-loops -floop-parallelize-all -fwhole-file -ftree-parallelize-loops=4 -floop-interchange -floop-strip-mine -floop-block")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -ffree-form -ffree-line-length-none -fexternal-blas -march=native -m64 -funroll-loops -fwhole-file -floop-interchange -floop-strip-mine -floop-block")

    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "-O0 --coverage")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "-O3")
    
    set (PROFILE_FLAGS "-p -g")


endif()

# Append user flags to debug and release variables
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${USER_COMPILE_FLAGS_DEBUG}")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${USER_COMPILE_FLAGS_RELEASE}")

# Append profile flags
if (PROFILE)
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${PROFILE_FLAGS}")
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${PROFILE_FLAGS}")
endif()




#--------------------------------------------------------------------------------------------------
#                                          Find HDF5
#--------------------------------------------------------------------------------------------------
# NOTE: export HDF5_ROOT path as shell environment variable such that 
#       $HDF5_ROOT contains the bin, lib, include directories for HDF5.
#       The find_package module looks first for $HDF5_ROOT.

find_package(HDF5 REQUIRED COMPONENTS Fortran HL Fortran_HL)
include_directories(${HDF5_INCLUDE_DIRS})


if (HDF5_FOUND)
    message(STATUS "HDF5_LIBRARIES: " ${HDF5_LIBRARIES})
endif()




#--------------------------------------------------------------------------------------------------
#                                         Find Lapack
#--------------------------------------------------------------------------------------------------
find_package(LAPACK REQUIRED)
message(STATUS "LAPACK_LIBRARIES: " ${LAPACK_LIBRARIES})





#--------------------------------------------------------------------------------------------------
#                                          Find MPI
#--------------------------------------------------------------------------------------------------

# add -DUSER_MPI=True to config to read environment variables
if (USER_MPI)
    message(STATUS "Reading environment variables for MPI_Fortran_COMPILER and MPI_Fortran_INCLUDE_PATH. Please be sure these are set as environment variables")
    set(MPI_Fortran_COMPILER $ENV{MPI_Fortran_COMPILER})
    set(MPI_Fortran_INCLUDE_PATH $ENV{MPI_Fortran_INCLUDE_PATH})
endif()


find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})









#---------------------------------------------------------------------------------------------------
#                                           Find METIS
#---------------------------------------------------------------------------------------------------

set(METIS_HINT_INCLUDE_DIR $ENV{METIS_INCLUDE_DIR})


find_path(METIS_INCLUDE_DIR metis.h HINTS ${METIS_HINT_INCLUDE_DIR})
find_library(METIS_LIBRARIES metis)

if(METIS_INCLUDE_DIR)
    if(METIS_LIBRARIES)
        set(METIS_FOUND "YES") 
    endif(METIS_LIBRARIES)
endif(METIS_INCLUDE_DIR)

message(STATUS "METIS INCLUDE DIR: " ${METIS_INCLUDE_DIR})
message(STATUS "METIS LIBRARIES: " ${METIS_LIBRARIES})
include_directories(${METIS_INCLUDE_DIR})









#--------------------------------------------------------------------------------------------------
#                                           Find Tecplot
#--------------------------------------------------------------------------------------------------
# Shouldn't need user-intervention here.
# First define TEC_ROOT path as cmake variable such that
# ${TEC_ROOT}/lib/libtecio.a  ${TEC_ROOT}/include/tecio.f90 are valid.
# ${TEC_ROOT} gets used by the find_package(Tecplot) module.

# Select TecIO library based on detected OS
if ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" ) 
    set(TEC_ROOT ${CMAKE_SOURCE_DIR}/tecio/mac)
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    set(TEC_ROOT ${CMAKE_SOURCE_DIR}/tecio/linux)
else()
    message(SEND_ERROR "Invalid operating system detected for TecIO")
endif()


find_package(Tecplot REQUIRED)
include_directories(${TECPLOT_INCLUDE_DIR})


message(STATUS "TECPLOT_LIBRARIES: " ${TECPLOT_LIBRARIES})
message(STATUS "TECPLOT_INCLUDE_DIR: " ${TECPLOT_INCLUDE_DIR})




#-------------------------------------------------------------------------------------------------
#                                   Add source tree and includes
#-------------------------------------------------------------------------------------------------
add_subdirectory(DNAD)                                                                        
add_subdirectory(src)

include_directories(${CMAKE_SOURCE_DIR}/src/tools)    # Need include here for macro-enhanced message handler messenger.h







#-------------------------------------------------------------------------------------------------
#                               Add core library and driver targets
#-------------------------------------------------------------------------------------------------
add_library(coredg STATIC ${LIB_SRCS})
#target_link_libraries(coredg ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${TECPLOT_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_link_libraries(coredg ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${TECPLOT_LIBRARIES} ${METIS_LIBRARIES})
set_target_properties(coredg PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")
                             #COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}"
                             #LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")

# Add driver executable linked to coredg library
add_executable(chidg ${EXE_SRCS})
target_link_libraries(chidg coredg ${HDF5_LIBRARIES} ${TECPLOT_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${METIS_LIBRARIES})
set_target_properties(chidg PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -lstdc++ -lpthread")
                            #COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}"
                            #LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS} -lstdc++ -lpthread")
                            # need to link stdc++ here for tecio
                            






#-------------------------------------------------------------------------------------------------
#                                       Install targets
#-------------------------------------------------------------------------------------------------
install(TARGETS chidg coredg
        RUNTIME DESTINATION /usr/local/chidg/bin
        LIBRARY DESTINATION /usr/local/chidg/lib
        ARCHIVE DESTINATION /usr/local/chidg/lib)














#---------------------------------------------------------------------------------------------------
#                               Testing:   enable with -DTESTING=True
#---------------------------------------------------------------------------------------------------
# To build tests, the unit-testing framework pFUnit must be compiled with MPI support.
# The environment variable PFUNIT must be set such that PFUNIT/gnu or PFUNIT/intel
# are valid.

if (TESTING)

    enable_testing()

    #
    # Set pFUnit installation based on compiler vendor
    #
    #if (${Fortran_VENDOR} MATCHES "INTEL")
    #    set(PF_ROOT $ENV{PFUNIT}/intel)
    #elseif (${Fortran_VENDOR} MATCHES "GNU")
    #    set(PF_ROOT $ENV{PFUNIT}/gnu)
    #endif()

    set(PF_ROOT $ENV{PFUNIT})

    #
    # Set location of pFUnit libs
    #
    set(PFU_LIB ${PF_ROOT}/lib)


    #
    # Add unit tests
    #
    add_subdirectory(test)


endif()






#-------------------------------------------------------------------------------------------------
#                                   Generate documentation
#-------------------------------------------------------------------------------------------------





# This turns off the deprecated $HDF5_DIR variable that shows up sometimes when running ccmake
# It probably depends on what version of FindHDF5.cmake is being used, but the correct variable
# to set is HDF5_ROOT
unset(HDF5_DIR CACHE)


#-------------------------------------------------------------------------------------------------
