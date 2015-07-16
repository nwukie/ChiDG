cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)
project(ucdg Fortran)

#------------------------------------------------------------------                                                                        
#              CMake macros for adding source files
#------------------------------------------------------------------                                                                        
include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)
include(cmake/macros/add_test_target.txt)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#------------------------------------------------------------------                                                                        
#                     Set compiler flags
#------------------------------------------------------------------                                                                        
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_COMPILE_FLAGS_DEBUG   "-O0 -debug -fpp -p -g -traceback -check all -check bounds -fp-stack-check -assume realloc_lhs")
#    set (CMAKE_Fortran_COMPILE_FLAGS_DEBUG   "-O3 -p -g")
#    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3")
    #set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3 -p -g -fpp")
    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3 -p -fpp -g -assume realloc_lhs")
    #set (CMAKE_Fortran_FLAGS_RELEASE "-O0 -g -traceback -check uninit -debug all")
    #set (CMAKE_Fortran_FLAGS_DEBUG "-p -g")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "-g -traceback -check all -check bounds")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

#if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # -fcray-pointer for tecio output
    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3 -fcray-pointer -ffree-form -ffree-line-length-none -cpp")
    #set (CMAKE_Fortran_FLAGS_RELEASE "-O0 -g -traceback -check uninit -debug all")
    #set (CMAKE_Fortran_FLAGS_DEBUG "-p -g")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")
endif()


#set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS_RELEASE})
#set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS_RELEASE})

set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS_DEBUG})
set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS_DEBUG})




# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


#------------------------------------------------------------------                                                                        
#                        Find HDF5
#------------------------------------------------------------------
# NOTE: export HDF5_ROOT path as shell environment variable such that 
#       $HDF5_ROOT/lib  $HDF5_ROOT/include are valid. The find_package
#       module looks first for $HDF5_ROOT
find_package(HDF5 COMPONENTS Fortran HL Fortran_HL)
include_directories(${HDF5_INCLUDE_DIRS})

message(STATUS "HDF5_LIBRARIES: " ${HDF5_LIBRARIES})



#------------------------------------------------------------------                                                                        
#                        Find Lapack
#------------------------------------------------------------------
find_package(LAPACK REQUIRED)



#------------------------------------------------------------------
#                        Find Tecplot
#------------------------------------------------------------------
# NOTE: export TEC_ROOT path as shell environment variables such that
#       $TEC_ROOT/lib  $TEC_ROOT/include are valid.
find_package(Tecplot REQUIRED)
include_directories(${TECPLOT_INCLUDE_DIR})

message(STATUS "TECPLOT_LIBRARIES: " ${TECPLOT_LIBRARIES})
message(STATUS "TECPLOT_INCLUDE_DIR: " ${TECPLOT_INCLUDE_DIR})




#------------------------------------------------------------------                                                                        
#                        Add source tree
#------------------------------------------------------------------
add_subdirectory(DNAD)                                                                        
add_subdirectory(src)




#------------------------------------------------------------------                                                                        
#              Add core library and driver targets
#------------------------------------------------------------------                                                                        
add_library(coredg STATIC ${LIB_SRCS})
target_link_libraries(coredg ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES} ${TECPLOT_LIBRARIES})
set_target_properties(coredg PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")

# Add driver program
add_executable(chidg ${EXE_SRCS})
target_link_libraries(chidg coredg ${HDF5_LIBRARIES} ${TECPLOT_LIBRARIES})
set_target_properties(chidg PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -lstdc++ -lgfortran")
                            # need to link lstdc++ here for tecio
                            # need to link to gfortran if compiling with ifort and 
                            # linking to atlas, which was compiled with gfortran



#------------------------------------------------------------------                                                                    
#                     Enable unit-testing                                   
#------------------------------------------------------------------                                                                        
set(UNITTEST TRUE CACHE BOOL "Enable building unit tests")
if (UNITTEST)
    enable_testing()
    set(PFUNIT /Users/nathanwukie/Documents/pbus/pFUnit/install/gnu/serial)
    set(PFU_LIB ${PFUNIT}/lib)

    # Add tests
    add_subdirectory(unittests)
endif()


#------------------------------------------------------------------                                                                        
#                    Generate documentation
#------------------------------------------------------------------                                                                        





# This turns off the deprecated $HDF5_DIR variable that shows up sometimes when running ccmake
# It probably depends on what version of FindHDF5.cmake is being used, but the correct variable
# to set is HDF5_ROOT
unset(HDF5_DIR CACHE)









