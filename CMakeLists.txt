project(chidg Fortran)
cmake_minimum_required(VERSION 3.7.1)

set(ChiDG_VERSION_MAJOR 0)
set(ChiDG_VERSION_MINOR 2)

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
include(cmake/macros/add_python_sources.txt)

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


# Assemble libraries in form for linking with f2py Python interface
# Accumulate low-level hdf5 libraries
foreach(LIBFULL ${HDF5_LIBRARIES})
    get_filename_component(LIB ${LIBFULL} NAME_WE)
    set(TRIMMED_LIBRARIES ${TRIMMED_LIBRARIES} ${LIB})
endforeach()

# Accumulate high-level hdf5 libraries
foreach(LIBFULL ${HDF5_HL_LIBRARIES})
    get_filename_component(LIB ${LIBFULL} NAME_WE)
    set(TRIMMED_LIBRARIES ${TRIMMED_LIBRARIES} ${LIB})
endforeach()

# Replace leading 'lib' with '-l' for each library 
# so we can use it as the link statement for f2py
foreach(LIB ${TRIMMED_LIBRARIES})
    string(REPLACE "lib" " -l" PROCESSED_NAME ${LIB})
    set(PYTHON_LINK_HDF5_LIBRARIES ${PYTHON_LINK_HDF5_LIBRARIES} ${PROCESSED_NAME})
endforeach()

# Get library path
foreach(LIB ${HDF5_LIBRARIES})
    get_filename_component(HDF5_LIBRARY_PATH ${LIB} DIRECTORY)
endforeach()


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



# Assemble libraries in form for linking with f2py Python interface
unset(TRIMMED_LIBRARIES)
# Trim path and extension
foreach(LIBFULL ${MPI_Fortran_LIBRARIES})
    get_filename_component(LIB ${LIBFULL} NAME_WE)
    set(TRIMMED_LIBRARIES ${TRIMMED_LIBRARIES} ${LIB})
endforeach()

# Replace leading 'lib' with '-l' for each library 
# so we can use it as the link statement for f2py
foreach(LIB ${TRIMMED_LIBRARIES})
    string(REPLACE "lib" " -l" PROCESSED_NAME ${LIB})
    set(PYTHON_LINK_MPI_LIBRARIES ${PYTHON_LINK_MPI_LIBRARIES} ${PROCESSED_NAME})
endforeach()

# Get library path for -L
foreach(LIB ${MPI_Fortran_LIBRARIES})
    get_filename_component(MPI_LIBRARY_PATH ${LIB} DIRECTORY)
endforeach()



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


    # Write test program to file.
    set(test_lapack_file ${CMAKE_CURRENT_BINARY_DIR}/test_lapack.f90)
    file(WRITE ${test_lapack_file}
        "program main\n"
        "    implicit none\n"
        "    external DGETRI\n"
        "    external DGETRF\n"
        "    real(8), dimension(:,:),    allocatable :: A(:,:), Ainv(:,:)\n"
        "    real(8), dimension(:),      allocatable :: work\n"
        "    integer,  dimension(:),     allocatable :: ipiv\n"
        "    integer :: n, info\n"
        "    allocate(A(20,20))\n"
        "    allocate(Ainv(20,20), work(20), ipiv(20))\n"
        "    Ainv = A\n"
        "    n = size(A,1)\n"
        "    call DGETRF(n, n, Ainv, n, ipiv, info)\n"
        "    call DGETRI(n, Ainv, n, ipiv, work, n, info)\n"
        "end program main\n")


    # Try compiling and running test program that uses BLAS/LAPACK routines.
    # We are checking here to make sure that the detected libraries actually
    # work on the current machine. 
    try_run(RUN_RESULTS COMPILE_RESULTS
            ${CMAKE_CURRENT_BINARY_DIR} ${test_lapack_file}
            LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
            OUTPUT_VARIABLE out_var)

    # Interpret RUN/COMPILE results as good/bad
    set(COMPILED_OKAY ${COMPILE_RESULTS})
    if (${RUN_RESULTS} EQUAL 0)
        set(RAN_OKAY TRUE)
    else()
        set(RAN_OKAY FALSE)
    endif()

    # Handle RUN/COMPILE results
    if (NOT ${COMPILED_OKAY})
        message(SEND_ERROR "We couldn't compile a test program using the BLAS/LAPACK \
                            Libraries that were detected. You can try and influence \
                            the BLAS library that is detected using the configuration\
                            parameter -DBLA_VENDOR. For example, to explicitly try\
                            to use the Apple Accelerate Framework you can pass\
                            -DBLA_VENDOR=Apple during configuration.")
    elseif(NOT ${RAN_OKAY} )
        message(SEND_ERROR "A test program failed when trying to use the BLAS/LAPACK \
                            Libraries that were detected. This could mean that they \
                            are not compatible with your machine, which is possible \
                            if you are using libraries from a network directory \
                            managed by another user. You can try and influence the \
                            BLAS library that is detected using the configuration \
                            parameter -DBLA_VENDOR. For example to explicitly try  \
                            to use the Apple Accelerate Framework you can pass \
                            -DBLA_VENDOR=Apple during configuration.")
    endif()


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
    install(FILES ${TEC_ROOT}/lib/libtecio.dylib DESTINATION /usr/local/chidg/lib)
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    set(TEC_ROOT ${CMAKE_SOURCE_DIR}/tecio/linux)
    install(FILES ${TEC_ROOT}/lib/libtecio.so DESTINATION /usr/local/chidg/lib)
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


include(cmake/macros/cmake_fortran_macros.txt)
get_fortran_type_size(type_size)
get_fortran_class_size(class_size)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/constants/mod_kinds.f90.in" "${CMAKE_CURRENT_BINARY_DIR}/mod_kinds.f90" @ONLY)
list(APPEND LIB_SRCS "${CMAKE_CURRENT_BINARY_DIR}/mod_kinds.f90")


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
                            


#
# Build python module
#
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(QUIET_STATUS --quiet)
endif()


# Copy 'pychidg' package to build directory so we can build the _chidg.so extension
# module inside of it and install it with a parameter in distutils package_data.
# This install behavior is set in the install script 'src/interfaces/setup.py.in'
# which is configured by cmake with some text substitutions for version numbers
# and directory locations.
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/src/interfaces/pychidg" 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


# PYTHON 1:
# Build python extension module using f2py to process opaque interfaces
set(PYTHON_INTERFACE_MODULE _chidg.so)
add_custom_command(OUTPUT ${PYTHON_INTERFACE_MODULE}
                   COMMAND f2py-f90wrap
                           -c ${PYTHON_SRCS}
                           -m _chidg
                           --f90exec=${CMAKE_Fortran_COMPILER}
                           --build-dir "${CMAKE_BINARY_DIR}/pychidg"
                           ${CMAKE_BINARY_DIR}/lib/libcoredg.a
                           ${METIS_LIBRARIES}
                           -L${TEC_ROOT}/lib
                           -ltecio
                           -L${HDF5_LIBRARY_PATH}
                           ${PYTHON_LINK_HDF5_LIBRARIES}
                           -L${MPI_LIBRARY_PATH}
                           ${PYTHON_LINK_MPI_LIBRARIES}
                           --link-lapack_opt
                           -I${CMAKE_BINARY_DIR}/lib
                           ${QUIET_STATUS}
                   DEPENDS ${PYTHON_SRCS} coredg
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pychidg
                   COMMENT "f2py: Building Python module _chidg")


# Define/configure setup.py 
set(SETUP_PY_IN  "${CMAKE_CURRENT_SOURCE_DIR}/src/interfaces/setup.py.in")
set(SETUP_PY     "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(DEPS         "${CMAKE_CURRENT_BINARY_DIR}/pychidg/__init__.py")
set(PYTHON_BUILD "${CMAKE_CURRENT_BINARY_DIR}/python-build-stamp")
configure_file(${SETUP_PY_IN} ${SETUP_PY})


# PYTHON 2:
# Build python module
find_program(PYTHON "python")
add_custom_command(OUTPUT ${PYTHON_BUILD}
                   COMMAND ${PYTHON} ${SETUP_PY} build
                   COMMAND ${CMAKE_COMMAND} -E touch ${PYTHON_BUILD}
                   DEPENDS ${DEPS})


# Add target 'make python'. Triggers sections PYTHON 1, PYTHON 2
add_custom_target(python DEPENDS ${PYTHON_INTERFACE_MODULE} ${PYTHON_BUILD})


# Define python install, if extension modules were built('make python')
# then execute package install script 'python setup.py install'
install(CODE 
        "if (EXISTS ${PYTHON_BUILD}) 
        execute_process(COMMAND ${PYTHON} ${SETUP_PY} install) 
        endif()"
        )


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
#   pFUnit from git submodule:
#   --------------------------
#   This is already set-up. One only needs to run 'make check' to build pFUnit and ChiDG tests.
#   This step requires internet access. If internet access is not available, one could provide
#   the pFUnit source as a .tar archive as detailed in the next section.
#   
#   pFUnit from .tar:
#   --------------------------
#   If PFUNIT_TAR is defined during configure 
#   cmake <other_stuff> -DPFUNIT_TAR=/home/user/pFUnit.tar <other_stuff>
#   Then this will bypass the git submodule approach for obtaining the pFUnit source.
#   This approach is useful if internet access is not available on the host machine.
#   
#   An archive of pFUnit can be obtained at: https://sourceforge.net/projects/pfunit
#
#----------------------------------------------------------------------------------------------
include(ExternalProject)


if (PFUNIT_TAR)
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory external
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                    )
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory pfunit-src
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external
                    )
    execute_process(COMMAND tar xzf ${PFUNIT_TAR} -C ${PROJECT_BINARY_DIR}/external/pfunit-src --strip-components=1
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                    )
    set(PFUNIT_DOWNLOAD_COMMAND "")
    set(PFUNIT_DOWNLOAD_DIR "")
    set(PFUNIT_SOURCE_DIR ${PROJECT_BINARY_DIR}/external/pfunit-src)

else()
    set(PFUNIT_DOWNLOAD_COMMAND git submodule update)
    set(PFUNIT_DOWNLOAD_DIR "${PROJECT_SOURCE_DIR}")
    set(PFUNIT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/pfunit")
endif()




set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DMPI=YES
    -DINSTALL_PATH=${PROJECT_BINARY_DIR}/external/pfunit   # keeps pfunit from complaining
    )



ExternalProject_Add(
    pfunit
    DOWNLOAD_COMMAND ${PFUNIT_DOWNLOAD_COMMAND}
    DOWNLOAD_DIR ${PFUNIT_DOWNLOAD_DIR}
    SOURCE_DIR ${PFUNIT_SOURCE_DIR}
    BINARY_DIR ${PROJECT_BINARY_DIR}/external/pfunit-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${ExternalProjectCMakeArgs}
    EXCLUDE_FROM_ALL 1
    )


if (NOT PFUNIT_TAR)
    add_dependencies(pfunit git_update)
endif()
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
