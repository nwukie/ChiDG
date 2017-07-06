project(chidg Fortran)
cmake_minimum_required(VERSION 3.7.1)


#-----------------------------------------------------------------------------------
#                           Check for bad CMake version
#-----------------------------------------------------------------------------------
if ( (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 6 AND CMAKE_PATCH_VERSION EQUAL 0) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 6 AND CMAKE_PATCH_VERSION EQUAL 1) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 7 AND CMAKE_PATCH_VERSION EQUAL 0) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 7 AND CMAKE_PATCH_VERSION EQUAL 2) )
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
#                               Set compiler flags
#----------------------------------------------------------------------------------- 
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -g \
                                     -traceback -check all -mcmodel=medium \
                                     -fp-stack-check -heap-arrays 1600")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -fpp -assume realloc_lhs -mcmodel=medium \
                                     -align array32byte -traceback")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

    set (PROFILE_FLAGS "-p")
    
    #set (CMAKE_AR xiar)        # for -ipo
    #set (CMAKE_LINKER xild)    # for -ipo
    
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)

    # Enforce gfortran minimum version requirement
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.0.0)
        message(FATAL_ERROR "Require at least gfortran 6.0.0.")
    endif()

    # Set compile/link flags
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp --coverage -ffree-form \
                                     -ffree-line-length-none -fbacktrace \
                                     -fstack-protector -fcheck=bounds \
                                     -p -g")    
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
find_package(HDF5 1.10 REQUIRED COMPONENTS Fortran HL Fortran_HL)
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

if (LAPACK_FOUND)
    message(STATUS "    LAPACK Libraries: ")
    foreach (L ${LAPACK_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
endif()





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
    foreach (L ${MPI_Fortran_INCLUDE_PATH})
        message(STATUS "        " ${L})
    endforeach()
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


#
# Initialize GIT_SHA1 variable with git commit hash
#
# Write git SHA1 to mod_version.f90. Text substitution in 'mod_version.f90.in' generates 'mod_version.f90'
#
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/constants/mod_version.f90.in" "${CMAKE_CURRENT_BINARY_DIR}/mod_version.f90" @ONLY)
list(APPEND LIB_SRCS "${CMAKE_CURRENT_BINARY_DIR}/mod_version.f90")



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
#
#                              Target to update git submodules
#
#----------------------------------------------------------------------------------------------
add_custom_target(
    git_update
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )






#----------------------------------------------------------------------------------------------
#
#   pFUnit external project for testing.
#
#   The following code is based on modifications of code by Randovan Bast. See file LICENSE, 
#   section BAST for license.
#   
#
#----------------------------------------------------------------------------------------------
include(ExternalProject)

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DMPI=YES
    -DINSTALL_PATH=${PROJECT_BINARY_DIR}/external/pfunit   # keeps pfunit from complaining
    )

ExternalProject_Add(
    pfunit
    DOWNLOAD_COMMAND git submodule update
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/pfunit
    BINARY_DIR ${PROJECT_BINARY_DIR}/external/pfunit-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${ExternalProjectCMakeArgs}
    EXCLUDE_FROM_ALL 1
    )


add_dependencies(pfunit git_update)
set(PF_ROOT ${PROJECT_BINARY_DIR}/external/pfunit)
set(PFU_LIB ${PF_ROOT}/lib)
set(MPIRUN_EXECUTABLE "mpirun")


#----------------------------------------------------------------------------------------------
#                               Testing:
#   To build tests:     make check
#   To run tests:       make test
#
#   To build tests, the unit-testing framework pFUnit must be compiled with MPI support.
#   The environment variable PFUNIT must be set such that PFUNIT/bin PFUNIT/lib
#   are valid.
#
#----------------------------------------------------------------------------------------------


#
# Add unit tests. Build with 'make check'. Run with 'make test'.
#
enable_testing()
add_custom_target(check)
add_subdirectory(test)



# This turns off the deprecated $HDF5_DIR variable that shows up sometimes when running ccmake
# It probably depends on what version of FindHDF5.cmake is being used, but the correct variable
# to set is HDF5_ROOT
unset(HDF5_DIR CACHE)


#----------------------------------------------------------------------------------------------
