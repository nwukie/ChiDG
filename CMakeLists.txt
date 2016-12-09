project(chidg Fortran)
cmake_minimum_required(VERSION 3.0)


#-----------------------------------------------------------------------------------
#                           Check for bad CMake version
#-----------------------------------------------------------------------------------
if ( (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 6 AND CMAKE_PATCH_VERSION EQUAL 0) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 6 AND CMAKE_PATCH_VERSION EQUAL 1) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 7 AND CMAKE_PATCH_VERSION EQUAL 0) )
    message(FATAL_ERROR "cmake versions 3.6.0, 3.6.1, 3.6.2, and 3.7.0 have regressions that \
                         break the build. Please use a different version of cmake.")
endif()



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
#                               Set default options
#-----------------------------------------------------------------------------------

# Set default options 

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build. Options are: Debug Release" FORCE)
endif()





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
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -p -g \
                                     -traceback -check all -check bounds \
                                     -fp-stack-check -standard-semantics \
                                     -assume nostd_mod_proc_name \
                                     -heap-arrays 1600")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -fpp -assume realloc_lhs")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

    set (PROFILE_FLAGS "")
    
    #set (CMAKE_AR xiar)        # for -ipo
    #set (CMAKE_LINKER xild)    # for -ipo
    
elseif (${Fortran_VENDOR} MATCHES "GNU")
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp --coverage -ffree-form \
                                     -ffree-line-length-none -fbacktrace \
                                     -fstack-protector -fcheck=bounds \
                                     -fmax-stack-var-size=1600 -p -g")    
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -ffree-form \
                                     -ffree-line-length-none -fexternal-blas \
                                     -funroll-loops")
                                     #-march=native -funroll-loops")

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




#----------------------------------------------------------------------------------------------
#                                          Find HDF5
#
#   NOTE: export HDF5_ROOT path as shell environment variable such that 
#         $HDF5_ROOT contains the bin, lib, include directories for HDF5.
#         The find_package module looks first for $HDF5_ROOT.
#
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "HDF5 Summary:")
find_package(HDF5 REQUIRED COMPONENTS Fortran HL Fortran_HL)
include_directories(${HDF5_INCLUDE_DIRS})


if (HDF5_FOUND)
    message(STATUS "    HDF5 Libraries: ")
    foreach (L ${HDF5_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
    message(STATUS "    HDF5 High-Level Libraries: ")
    foreach (L ${HDF5_HL_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
endif()




#----------------------------------------------------------------------------------------------
#                                         Find Lapack
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "BLAS/LAPACK Summary:")
find_package(LAPACK REQUIRED)
message(STATUS "    LAPACK Libraries: ")
message(STATUS "        " ${LAPACK_LIBRARIES})





#----------------------------------------------------------------------------------------------
#                                          Find MPI
#
#   NOTE: It is important to get MPI libraries + modules were built with a similar compiler
#         
#   
#   If MPI is not detected by default, or a specific implementation needs forced, set
#   the environment variable: MPI_Fortran_COMPILER
#
#       ex.(bash): export MPI_Fortran_COMPILER=/usr/local/bin/mpifort
#
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "MPI Summary:")
set(MPI_Fortran_COMPILER $ENV{MPI_Fortran_COMPILER})



find_package(MPI REQUIRED COMPONENTS Fortran)
include_directories(${MPI_Fortran_INCLUDE_PATH})



if (MPI_Fortran_FOUND)
    message(STATUS "    MPI Fortran Libraries: ")
    foreach (L ${MPI_Fortran_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
    message(STATUS "    MPI Include Path: ")
    message(STATUS "        " ${MPI_Fortran_INCLUDE_PATH})
endif()







#----------------------------------------------------------------------------------------------
#                                           Find METIS
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "METIS Summary: ")
find_package(METIS)
include_directories(${METIS_INCLUDE_DIR})


if(NOT ${METIS_FOUND})
    message(SEND_ERROR "METIS library not found")
else()
    message(STATUS "    METIS Libraries:")
    message(STATUS "        " ${METIS_LIBRARIES})
endif()









#----------------------------------------------------------------------------------------------
#                                           Find Tecplot
# Shouldn't need user-intervention here.
#
# TEC_ROOT is defined as a path such that 
#       ${TEC_ROOT}/lib/libtecio.a  
#       ${TEC_ROOT}/include/tecio.f90 
# are valid.
#
# ${TEC_ROOT} gets used by the find_package(Tecplot) module.
#
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "Tecplot Summary: ")

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


message(STATUS "    Tecplot Libraries: ")
message(STATUS "        " ${TECPLOT_LIBRARIES})


#----------------------------------------------------------------------------------------------
#                                   Add source tree and includes
#----------------------------------------------------------------------------------------------
add_subdirectory(DNAD)                                                                        
add_subdirectory(src)

# Need include here for macro-enhanced message handler messenger.h
include_directories(${CMAKE_SOURCE_DIR}/src/tools)







#----------------------------------------------------------------------------------------------
#                                         Build library
#----------------------------------------------------------------------------------------------
message(STATUS " ")
#add_library(coredg SHARED ${LIB_SRCS})
add_library(coredg STATIC ${LIB_SRCS})
target_link_libraries(coredg 
                      ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} 
                      ${LAPACK_LIBRARIES}
                      ${BLAS_LIBRARIES} 
                      ${TECPLOT_LIBRARIES} 
                      ${MPI_Fortran_LIBRARIES}
                      ${METIS_LIBRARIES}
                      )
set_target_properties(coredg PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")





#----------------------------------------------------------------------------------------------
#                                       Build Interfaces
#----------------------------------------------------------------------------------------------
# Add driver executable linked to coredg library
add_executable(chidg ${EXE_SRCS})
target_link_libraries(chidg coredg)
set_target_properties(chidg PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")
                            




#----------------------------------------------------------------------------------------------
#                                       Install targets
#----------------------------------------------------------------------------------------------
install(TARGETS chidg coredg
        RUNTIME DESTINATION /usr/local/chidg/bin
        LIBRARY DESTINATION /usr/local/chidg/lib
        ARCHIVE DESTINATION /usr/local/chidg/lib)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION /usr/local/chidg/include)












#----------------------------------------------------------------------------------------------
#                               Testing:
#   To build tests:     make check
#   To run tests:       make test
#
#----------------------------------------------------------------------------------------------
# To build tests, the unit-testing framework pFUnit must be compiled with MPI support.
# The environment variable PFUNIT must be set such that PFUNIT/gnu or PFUNIT/intel
# are valid.


file(REMOVE ${CMAKE_BINARY_DIR}/bin/unittest.sh)
file(REMOVE ${CMAKE_BINARY_DIR}/bin/regtest.sh)


enable_testing()

add_custom_target(check)

#
# Get paths for pFUnit and MPIRUN
#
set(PF_ROOT $ENV{PFUNIT})
set(MPIRUN_EXECUTABLE $ENV{MPIRUN_EXECUTABLE})


#
# Set location of pFUnit libs
#
set(PFU_LIB ${PF_ROOT}/lib)


#
# Add unit tests
#
add_subdirectory(test)









# This turns off the deprecated $HDF5_DIR variable that shows up sometimes when running ccmake
# It probably depends on what version of FindHDF5.cmake is being used, but the correct variable
# to set is HDF5_ROOT
unset(HDF5_DIR CACHE)


#----------------------------------------------------------------------------------------------
