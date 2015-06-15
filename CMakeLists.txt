cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)
project(ucdg Fortran)

#------------------------------------------------------------------                                                                        
#              CMake macros for adding source files
#------------------------------------------------------------------                                                                        
include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)
include(cmake/macros/add_test_target.txt)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#------------------------------------------------------------------                                                                        
#                     Set compiler flags
#------------------------------------------------------------------                                                                        
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_COMPILE_FLAGS_DEBUG   "-O0 -debug -p -g -traceback -check all -check bounds -fp-stack-check")
#    set (CMAKE_Fortran_COMPILE_FLAGS_DEBUG   "-O3 -p -g")
#    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3")
    set (CMAKE_Fortran_COMPILE_FLAGS_RELEASE "-O3 -p -g -fpp")
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


set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS_RELEASE})
set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS_RELEASE})

#set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS_DEBUG})
#set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS_DEBUG})




# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


#------------------------------------------------------------------                                                                        
#                        Find HDF5
#------------------------------------------------------------------
set(HDF5_LIBRARY_DIR "~/Documents/pbus/hdf/hdf5-1.8.15/build/lib")
set(HDF5_INCLUDE_PATH "~/Documents/pbus/hdf/hdf5-1.8.15/build/include")
include_directories(${HDF5_INCLUDE_PATH})
find_library(HDF5_LIB NAMES libhdf5.a PATHS ${HDF5_LIBRARY_DIR} NO_DEFAULT_PATH)
find_library(HDF5_HL NAMES libhdf5_hl.a PATHS ${HDF5_LIBRARY_DIR} NO_DEFAULT_PATH)
find_library(HDF5_FHL NAMES libhdf5hl_fortran.a PATHS ${HDF5_LIBRARY_DIR} NO_DEFAULT_PATH)
find_library(HDF5_F NAMES libhdf5_fortran.a PATHS ${HDF5_LIBRARY_DIR} NO_DEFAULT_PATH)

message(STATUS "HDF5_LIBRARY: " ${HDF5_LIB})
message(STATUS "HDF5_LIBRARY: " ${HDF5_HL})
message(STATUS "HDF5_LIBRARY: " ${HDF5_FHL})
message(STATUS "HDF5_LIBRARY: " ${HDF5_F})

find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    #target_link_libraries(MyProg ${ZLIB_LIBRARIES})
endif()



#------------------------------------------------------------------                                                                        
#                        Find CGNS
#------------------------------------------------------------------
#set(CGNS_LIBRARY_DIR "~/Documents/pbus/cgns/cgnslib_3.2.1/build_hdf5/lib")
#set(CGNS_INCLUDE_PATH "~/Documents/pbus/cgns/cgnslib_3.2.1/build_hdf5/include")
#find_path(CGNS_INCLUDE_DIR cgnslib.h PATHS ${CGNS_INCLUDE_PATH} NO_DEFAULT_PATH)
#find_library(CGNS_LIBRARIES NAMES libcgns.a PATHS ${CGNS_LIBRARY_DIR} NO_DEFAULT_PATH)
#message(STATUS "CGNS_LIBRARIES: " ${CGNS_LIBRARIES})
#set(CGNS_LIBRARIES ${CGNS_LIBRARIES} ${HDF5_LIBRARIES})


#------------------------------------------------------------------                                                                        
#                        Add source tree
#------------------------------------------------------------------                                                                        
add_subdirectory(src)




#------------------------------------------------------------------                                                                        
#              Add core library and driver targets
#------------------------------------------------------------------                                                                        
add_library(coredg STATIC ${LIB_SRCS})
target_link_libraries(coredg ${HDF5_LIB} ${HDF5_HL} ${HDF5_FHL} ${HDF5_F} ${ZLIB_LIBRARIES})
set_target_properties(coredg PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")

# Add driver program
add_executable(flexdg ${EXE_SRCS})
target_link_libraries(flexdg coredg ${HDF5_LIB} ${HDF5_HL} ${HDF5_FHL} ${HDF5_F})
set_target_properties(flexdg PROPERTIES 
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
    set(PFUNIT /Users/nathanwukie/Documents/pbus/pFUnit/install/gnu/serial)
    set(PFU_LIB ${PFUNIT}/lib)


    # Add tests
    add_subdirectory(unittests)
endif()


#------------------------------------------------------------------                                                                        
#                    Generate documentation
#------------------------------------------------------------------                                                                        










