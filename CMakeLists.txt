cmake_minimum_required(VERSION 3.1)
enable_language(Fortran)
project(ucdg Fortran)

#------------------------------------------------------------------                                                                        
#              CMake macros for adding source files
#------------------------------------------------------------------                                                                        
include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)
include(cmake/macros/add_test_target.txt)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)



#------------------------------------------------------------------ 
#                    Set 'Debug' and 'Release' options
#------------------------------------------------------------------
set(UNITTEST FALSE CACHE BOOL "Enable building unit tests")
if (${CMAKE_BUILD_TYPE} MATCHES "Debug") 
    set(UNITTEST TRUE)
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(UNITTEST FALSE)
endif()





#------------------------------------------------------------------
#                    Get fortran compiler vendor
#------------------------------------------------------------------
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (Fortran_VENDOR "INTEL")
elseif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set (Fortran_VENDOR "GNU")
else ()
    message(SEND_ERROR "Fortran compiler not recognized")
endif()


#------------------------------------------------------------------                                                                        
#                     Set compiler flags
#------------------------------------------------------------------                                                                        
if (${Fortran_VENDOR} MATCHES "INTEL")
#    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -p -g -traceback -check all -check noarg_temp_created -check bounds -fp-stack-check")
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -p -g -traceback -check all -check bounds -fp-stack-check -standard-semantics -assume nostd_mod_proc_name")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -fpp -assume realloc_lhs")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")
    
    #set (CMAKE_AR xiar)        # for -ipo
    #set (CMAKE_LINKER xild)    # for -ipo
    
elseif (${Fortran_VENDOR} MATCHES "GNU")
    #set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp -ffree-form -ffree-line-length-none -fcheck=all -fbacktrace -fcheck=bounds -p -g")    # Note: -fcheck=all has issues with passing null pointer in tecio_interface
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp -ffree-form -ffree-line-length-none -fbacktrace -fcheck=bounds -p -g")    # Note: -fcheck=all has issues with passing null pointer in tecio_interface
#    set (USER_COMPILE_FLAGS_DEBUG "-O3 -cpp -ffree-form -ffree-line-length-none")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -ffree-form -ffree-line-length-none")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")
    
endif()

# Append user flags to debug and release variables
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${USER_COMPILE_FLAGS_DEBUG}")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${USER_COMPILE_FLAGS_RELEASE}")






#------------------------------------------------------------------                                                                        
#                        Find HDF5
#------------------------------------------------------------------
# NOTE: export HDF5_ROOT path as shell environment variable such that 
#       $HDF5_ROOT/gnu or $HDF5_ROOT/intel are valid. Inside the 'gnu' or 'intel'
#       folder one should find lib and include directories. The find_package
#       module looks first for $HDF5_ROOT. Here we append the gnu or intel folder to 
#       the environment variable based on what compiler we are using.
if (${Fortran_VENDOR} MATCHES "INTEL")
    message(STATUS "$ENV{HDF5_ROOT}")
    set(ENV{HDF5_ROOT} $ENV{HDF5_ROOT}/intel)
    message(STATUS "$ENV{HDF5_ROOT}")
    find_package(HDF5 REQUIRED COMPONENTS Fortran HL Fortran_HL)
#    find_package(HDF5 REQUIRED COMPONENTS Fortran HL Fortran_HL HINTS intel)
    include_directories(${HDF5_INCLUDE_DIRS})
elseif (${Fortran_VENDOR} MATCHES "GNU")
        message(STATUS "$ENV{HDF5_ROOT}")
    set(ENV{HDF5_ROOT} $ENV{HDF5_ROOT}/gnu)
    message(STATUS "$ENV{HDF5_ROOT}")
    find_package(HDF5 REQUIRED COMPONENTS Fortran HL Fortran_HL)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

if (HDF5_FOUND)
    message(STATUS "HDF5_LIBRARIES: " ${HDF5_LIBRARIES})
endif()



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
#                        Add source tree and includes
#------------------------------------------------------------------
add_subdirectory(DNAD)                                                                        
add_subdirectory(src)
include_directories(src/tools)    # Need include here for macro-enhanced message handler message.h





#------------------------------------------------------------------                                                                        
#              Add core library and driver targets
#------------------------------------------------------------------                                                                        
add_library(coredg STATIC ${LIB_SRCS})
target_link_libraries(coredg ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES} ${TECPLOT_LIBRARIES})
set_target_properties(coredg PROPERTIES 
                             COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                             LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")

# Add driver program
add_executable(chidg ${EXE_SRCS})
target_link_libraries(chidg coredg ${HDF5_LIBRARIES} ${TECPLOT_LIBRARIES})
set_target_properties(chidg PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -lstdc++ -lpthread")
                            # need to link lstdc++ here for tecio
                            

#------------------------------------------------------------------                                                                    
#                     Unit-testing                                   
#------------------------------------------------------------------                                                                        
if (UNITTEST)
    enable_testing()
    if (${Fortran_VENDOR} MATCHES "INTEL")
        set(PF_ROOT $ENV{PFUNIT}/intel)
    elseif (${Fortran_VENDOR} MATCHES "GNU")
        set(PF_ROOT $ENV{PFUNIT}/gnu)
    endif()
    set(PFU_LIB ${PF_ROOT}/lib)

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









