project(chidg Fortran C CXX)
cmake_minimum_required(VERSION 3.7.1)
include(ExternalProject)

set(ChiDG_VERSION_MAJOR 0)
set(ChiDG_VERSION_MINOR 2)

#-----------------------------------------------------------------------------------
#                           Check for bad CMake version
#-----------------------------------------------------------------------------------
if ( (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 6 AND CMAKE_PATCH_VERSION EQUAL 0) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 6 AND CMAKE_PATCH_VERSION EQUAL 1) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 7 AND CMAKE_PATCH_VERSION EQUAL 0) OR
     (CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION EQUAL 7 AND CMAKE_PATCH_VERSION EQUAL 2) )
    message(FATAL_ERROR "cmake versions 3.6.0, 3.6.1, 3.6.2, 3.7.0, and 3.7.2 have regressions that \
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
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.3.0)
        message(FATAL_ERROR "Require at least gfortran 6.3.0.")
    endif()

    # Set compile/link flags
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp -ffree-form \
                                     -ffree-line-length-none -fbacktrace \
                                     -fstack-protector -fcheck=bounds \
                                     -g")    
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -ffree-form \
                                     -ffree-line-length-none -fexternal-blas \
                                     -funroll-loops")
                                     #-march=native -funroll-loops")

    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "-O0")
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
set(MPI_C_COMPILER $ENV{MPI_C_COMPILER})
set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})



#find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(MPI REQUIRED COMPONENTS CXX Fortran)
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${MPI_CXX_INCLUDE_PATH})

message(STATUS "        " ${MPI_Fortran_LIBRARIES})
message(STATUS "        " ${MPI_CXX_LIBRARIES})

if (MPI_Fortran_FOUND)
    message(STATUS "    MPI Fortran Libraries: ")
    foreach (L ${MPI_Fortran_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
    message(STATUS "    MPI Fortran Include Path: ")
    foreach (L ${MPI_Fortran_INCLUDE_PATH})
        message(STATUS "        " ${L})
    endforeach()
endif()

if (MPI_CXX_FOUND)
    message(STATUS "    MPI CXX Libraries: ")
    foreach (L ${MPI_CXX_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
    message(STATUS "    MPI CXX Include Path: ")
    foreach (L ${MPI_CXX_INCLUDE_PATH})
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
foreach(LIBFULL ${MPI_CXX_LIBRARIES})
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
    set(LAPACK_DIRECTORIES)
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
    #
    # NOTE: 
    # We only try this when the project is first generated. For some reason, 
    # if CMake needs to regenerate, it seems to remove the rpath in the compiled
    # program, which then fails because it can't load the blas/lapack libraries.
    # This was a problem when compiling with MKL, but not Apple Accelerate.
    #
    set(LAPACK_TESTED FALSE CACHE BOOL "")
    mark_as_advanced(LAPACK_TESTED)
    if (NOT ${LAPACK_TESTED})
        try_run(RUN_RESULTS COMPILE_RESULTS
                ${CMAKE_CURRENT_BINARY_DIR} ${test_lapack_file}
                LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                COMPILE_OUTPUT_VARIABLE compile_out
                RUN_OUTPUT_VARIABLE run_out)
    endif()


    # Interpret RUN/COMPILE results as good/bad
    if (${RUN_RESULTS} EQUAL 0)
        set(RAN_OKAY TRUE)
    else()
        set(RAN_OKAY FALSE)
    endif()

    # Handle RUN/COMPILE results
    if (NOT ${COMPILED_OKAY})
        message(STATUS ${compile_out} )
        message(SEND_ERROR "We couldn't compile a test program using the BLAS/LAPACK \
                            Libraries that were detected. You can try and influence \
                            the BLAS library that is detected using the configuration\
                            parameter -DBLA_VENDOR. For example, to explicitly try\
                            to use the Apple Accelerate Framework you can pass\
                            -DBLA_VENDOR=Apple during configuration.")
    elseif(NOT ${RAN_OKAY} )
        message(STATUS ${run_out} )
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

    # If we made it through, signal that LAPACK has been tested so we don't need 
    # to test it again, since tested in regeneration seems to fail due to rpath issues.
    set(LAPACK_TESTED TRUE CACHE BOOL "" FORCE)

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
#   TecIO MPI uses some of the header-only Boost libraries. These have
#   been extracted from the Boost source into a Boost subset that includes
#   only the necessary Boost files that TecIO MPI depends on. These are
#   archived in tecio/tecio_boost.tar.gz.
#
#   The TecIO MPI source is distributed here in tecio/teciompisrc. This is build as 
#   a separate project and statically linked with the ChiDG executable.
#
#   Step 1: Unpack tecio_boost.tar
#   Step 2: Configure and build teciompi
#   Step 3: Set libraries variable, include directory to later link with ChiDG.
#
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "Tecplot Summary: ")


# Unpack Boost headers to support building teciompi
set(BOOST_UNPACKED FALSE CACHE BOOL "")
mark_as_advanced(BOOST_UNPACKED)
if (NOT ${BOOST_UNPACKED})
    message(STATUS "    Unpacking tecio_boost.tar.gz ...")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/tecio/tecio_boost.tar.gz
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
    message(STATUS "    Done: ${CMAKE_BINARY_DIR}/tecio_boost")
    set(BOOST_UNPACKED TRUE CACHE BOOL "" FORCE)
else()
    message(STATUS "    Boost already unpacked: ${CMAKE_BINARY_DIR}/tecio_boost")
endif()


set(TecIOMPICMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/teciompi
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DBOOST_ROOT=${CMAKE_BINARY_DIR}/tecio_boost
    )


ExternalProject_Add(
    teciompi
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/tecio/teciompisrc
    BINARY_DIR ${PROJECT_BINARY_DIR}/teciompi/teciompi-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/teciompi/teciompi-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/teciompi/teciompi-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/bin
    INSTALL_COMMAND cmake -E echo "Skipping install step."
    CMAKE_ARGS ${TecIOMPICMakeArgs}
    EXCLUDE_FROM_ALL 1
    )


set(TECPLOT_LIBRARIES ${PROJECT_BINARY_DIR}/teciompi/teciompi-build/libteciompi.a "stdc++")
include_directories(${PROJECT_SOURCE_DIR}/tecio/teciompisrc)





#----------------------------------------------------------------------------------------------
#
#   TIOGA external project for hole-cutting.
#
#   The following code is based on modifications of code by Randovan Bast. See file LICENSE, 
#   section BAST for license.
#
#   TIOGA from git submodule:
#   --------------------------
#   This is already set-up. This step requires internet access to clone the TIOGA project 
#   from git. If internet access is not available, one can provide the TIOGA source as 
#   a .tar archive as detailed in the next section.
#   
#   TIOGA from .tar:
#   --------------------------
#   If TIOGA_TAR is defined during configure:
#       cmake <other_stuff> -DTIOGA_TAR=/example/path/tioga.tar <other_stuff>
#   Then this will bypass the git submodule approach for obtaining the TIOGA source.
#   This approach is useful if internet access is not available on the host machine.
#   
#   An archive of TIOGA can be obtained at: 
#       https://github.com/jsitaraman/tioga
#
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "TIOGA Summary: ")


if (TIOGA_TAR)
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory external
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                    )
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory tioga-src
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external
                    )
    execute_process(COMMAND tar xzf ${TIOGA_TAR} -C ${PROJECT_BINARY_DIR}/external/tioga-src --strip-components=1
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                    )
    set(TIOGA_DOWNLOAD_COMMAND "")
    set(TIOGA_DOWNLOAD_DIR "")
    set(TIOGA_SOURCE_DIR ${PROJECT_BINARY_DIR}/external/tioga-src)

else()
    set(TIOGA_DOWNLOAD_COMMAND git submodule update)
    set(TIOGA_DOWNLOAD_DIR "${PROJECT_SOURCE_DIR}")
    set(TIOGA_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/tioga")
endif()




set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/tioga
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DTIOGA_HAS_NODEGID=OFF
    -DBUILD_TIOGA_EXE=OFF
    -DBUILD_GRIDGEN_EXE=OFF
    )



ExternalProject_Add(
    tioga
    DOWNLOAD_COMMAND ${TIOGA_DOWNLOAD_COMMAND}
    DOWNLOAD_DIR ${TIOGA_DOWNLOAD_DIR}
    SOURCE_DIR ${TIOGA_SOURCE_DIR}
    BINARY_DIR ${PROJECT_BINARY_DIR}/external/tioga-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/tioga-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/tioga-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${ExternalProjectCMakeArgs}
    )

if (NOT TIOGA_TAR)
    add_dependencies(tioga git_update)
endif()

set(TIOGA_LIBRARIES ${PROJECT_BINARY_DIR}/external/tioga/lib/libtioga.dylib)
set(TIOGA_INCLUDE_DIR ${PROJECT_BINARY_DIR}/external/tioga/include)
include_directories(${TIOGA_INCLUDE_DIR})

message(STATUS "    TIOGA Libraries: "     ${TIOGA_LIBRARIES})


#----------------------------------------------------------------------------------------------
#                                           Find PETSc
#----------------------------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "PETSc Summary:")
set(PETSC_EXECUTABLE_RUNS YES) # Override variable, since the module for petsc does not yet treat Fortran properly.
find_package (PETSc 3.11.3 REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})


if (PETSC_FOUND)
    message(STATUS "    PETSC Libraries: ")
    foreach (L ${PETSC_LIBRARIES})
        message(STATUS "        " ${L})
    endforeach()
    message(STATUS "    PETSC Compiler: ")
    message(STATUS "        " ${PETSC_COMPILER})



    # Write test program to file.
    set(test_petsc_file ${CMAKE_CURRENT_BINARY_DIR}/test_petsc.f90)

    # Try compiling and running test
    # petsc/src/vec/vec/examples/tutorials/ex11f.F
    file(WRITE ${test_petsc_file}
        "program main\n"
        "#include <petsc/finclude/petscvec.h>\n"
        "   use petscvec\n"
        "   implicit none\n"
        "   Vec               x\n"
        "   PetscReal         norm\n"
        "   PetscBool         flg\n"
        "   PetscMPIInt       rank\n"
        "   PetscInt          n,bs,comp\n"
        "   PetscErrorCode    ierr\n"
        "   PetscScalar       one\n"
        "   call PetscInitialize(PETSC_NULL_CHARACTER,ierr)\n"
        "   print*, ierr\n"
        "   if (ierr .ne. 0) stop\n"
        "   call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)\n"
        "   n   = 20\n"
        "   one = 1.0\n"
        "   call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, '-n',n,flg,ierr)\n"
        "   call VecCreate(PETSC_COMM_WORLD,x,ierr)\n"
        "   call VecSetSizes(x,PETSC_DECIDE,n,ierr)\n"
        "   bs = 2\n"
        "   call VecSetBlockSize(x,bs,ierr)\n"
        "   call VecSetFromOptions(x,ierr)\n"
        "   call VecSet(x,one,ierr)\n"
        "   call VecNorm(x,NORM_2,norm,ierr)\n"
        "   comp = 0\n"
        "   call VecStrideNorm(x,comp,NORM_2,norm,ierr)\n"
        "   comp = 1\n"
        "   call VecStrideNorm(x,comp,NORM_2,norm,ierr)\n"
        "   call VecStrideNorm(x,comp,NORM_1,norm,ierr)\n"
        "   call VecStrideNorm(x,comp,NORM_INFINITY,norm,ierr)\n"
        "   call VecDestroy(x,ierr)\n"
        "   call PetscFinalize(ierr)\n"
        "end program main\n")


    # Try compiling and running test program that uses PETSC routines.
    # We are checking here to make sure that the detected libraries actually
    # work on the current machine. 
    set(PETSC_TESTED FALSE CACHE BOOL "")
    mark_as_advanced(PETSC_TESTED)
    if (NOT ${PETSC_TESTED})
        try_run(RUN_RESULTS COMPILE_RESULTS
                ${CMAKE_CURRENT_BINARY_DIR} ${test_petsc_file}
                CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${PETSC_INCLUDES}"
                COMPILE_DEFINITIONS ${CMAKE_Fortran_FLAGS_RELEASE}
                LINK_LIBRARIES ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_Fortran_LIBRARIES}
                COMPILE_OUTPUT_VARIABLE compile_out
                RUN_OUTPUT_VARIABLE run_out)
    endif()

    # Interpret RUN/COMPILE results as good/bad
    set(COMPILED_OKAY ${COMPILE_RESULTS})
    if (${RUN_RESULTS} EQUAL 0)
        set(RAN_OKAY TRUE)
    else()
        set(RAN_OKAY FALSE)
    endif()

    # Handle RUN/COMPILE results
    if (NOT ${COMPILED_OKAY})
        message(STATUS ${compile_out} )
        message(SEND_ERROR "We couldn't compile a test program using the PETSC \
                            Libraries that were detected.")
    elseif(NOT ${RAN_OKAY} )
        message(STATUS ${run_out} )
        message(SEND_ERROR "A test program compiled but failed when trying to use the PETSC \
                            Libraries that were detected.")
    endif()


    # If we made it through, signal that PETSc has been tested so we don't need 
    # to test it again, since tested in regeneration seems to fail due to rpath issues.
    set(PETSC_TESTED TRUE CACHE BOOL "" FORCE)


endif()





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
add_library(chidg_library STATIC ${LIB_SRCS})
add_dependencies(chidg_library teciompi tioga)
target_link_libraries(chidg_library
                      ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} 
                      ${LAPACK_LIBRARIES}
                      ${BLAS_LIBRARIES} 
                      ${TECPLOT_LIBRARIES} 
                      ${MPI_Fortran_LIBRARIES}
                      ${METIS_LIBRARIES}
                      ${PETSC_LIBRARIES}
                      ${TIOGA_LIBRARIES}
                      )
set_target_properties(chidg_library PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}"
                             OUTPUT_NAME   "chidg")






#----------------------------------------------------------------------------------------------
#                                       Build Interfaces
#----------------------------------------------------------------------------------------------
# Add driver executable linked to chidg library
add_executable(chidg_binary ${EXE_SRCS})
target_link_libraries(chidg_binary chidg_library)
set_target_properties(chidg_binary PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}"
                            OUTPUT_NAME   "chidg")
                            


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
                           ${CMAKE_BINARY_DIR}/lib/libchidg.a
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
                   DEPENDS ${PYTHON_SRCS} chidg
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
install(TARGETS chidg_binary chidg_library
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/chidg)



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
