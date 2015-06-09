cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)
project(ucdg Fortran)


get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#=============================================
# Set compiler flags
#=============================================
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_COMPILE_FLAGS_DEBUG   "-O0 -debug -p -g -traceback -check all -check bounds -fp-stack-check")
#    set (CMAKE_Fortran_COMPILE_FLAGS_DEBUG   "-O3 -p -g")
#    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3")
    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3 -p -g")
    #set (CMAKE_Fortran_FLAGS_RELEASE "-O0 -g -traceback -check uninit -debug all")
    #set (CMAKE_Fortran_FLAGS_DEBUG "-p -g")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "-g -traceback -check all -check bounds")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

#if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # -fcray-pointer for tecio output
    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3 -fcray-pointer -ffree-form -ffree-line-length-none")
    #set (CMAKE_Fortran_FLAGS_RELEASE "-O0 -g -traceback -check uninit -debug all")
    #set (CMAKE_Fortran_FLAGS_DEBUG "-p -g")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")
endif()


set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS_RELEASE})
set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS_RELEASE})

#set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS_DEBUG})
#set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS_DEBUG})


#============================================
# MACROS for adding sources
#============================================
include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)


# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


#============================================
# Add source directory
#============================================
add_subdirectory(src)


# Add core library
add_library(coredg STATIC ${LIB_SRCS})
set_target_properties(coredg PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")

# Add driver program
add_executable(flexdg ${EXE_SRCS})
target_link_libraries(flexdg coredg)
set_target_properties(flexdg PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -lstdc++ -lgfortran")
                            # need to link lstdc++ here for tecio
                            # need to link to gfortran if compiling with ifort and 
                            # linking to atlas, which was compiled with gfortran



#=======================================================
#    UNIT TESTS
#=======================================================
set(UNITTEST FALSE CACHE BOOL "Enable building unit tests")
if (UNITTEST)
    set(PFUNIT /Users/nathanwukie/Documents/pbus/pFUnit/install/gnu/serial)
    set(PFU_LIB ${PFUNIT}/lib)


    # Add tests
    add_subdirectory(tests)
endif()


