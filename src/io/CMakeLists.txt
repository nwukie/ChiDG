add_lib_sources(mod_io.f90
                mod_hdfio.f90
                mod_hdf_utilities.f90
                mod_tecio_interface.f90
                mod_tecio.f90
                type_file_properties.f90
                mod_file_utilities.f90
                )












#----------------------------------------------------------------------------------------------------------
#                                   Set compiler flags for 'PLOT3DtoHDF5
#----------------------------------------------------------------------------------------------------------
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (COMPILE_FLAGS_plot3dtohdf5_DEBUG   "-O0 -debug -p -g -traceback -check all -check bounds -fp-stack-check -assume realloc_lhs")
    set (COMPILE_FLAGS_plot3dtohdf5_RELEASE "-O3 -assume realloc_lhs")
    set (LINK_FLAGS_plot3dtohdf5_DEBUG "-g -traceback -check all -check bounds")
    set (LINK_FLAGS_plot3dtohdf5_RELEASE "")

elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

    set (COMPILE_FLAGS_plot3dtohdf5_DEBUG "-O0 -fcray-pointer -ffree-form -ffree-line-length-none -cpp")
    set (COMPILE_FLAGS_plot3dtohdf5_RELEASE "-O3 -fcray-pointer -ffree-form -ffree-line-length-none -cpp")
    set (LINK_FLAGS_plot3dtohdf5_DEBUG "")
    set (LINK_FLAGS_plot3dtohdf5_RELEASE "")
endif()



add_executable(PLOT3DtoHDF5 plot3dtohdf5.f90)
target_link_libraries(PLOT3DtoHDF5 coredg ${HDF5_LIBRARIES})
set_target_properties(PLOT3DtoHDF5 PROPERTIES 
                            COMPILE_FLAGS "${COMPILE_FLAGS_plot3dtohdf5_RELEASE}"
                            LINK_FLAGS    "${LINK_FLAGS_plot3dtohdf5_RELEASE}")





#-----------------------------------------------------------------------------------------------------------
#                                   Set compiler flags for 'HDF5toTEC'
#-----------------------------------------------------------------------------------------------------------
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (COMPILE_FLAGS_HDF5toTEC_DEBUG   "-O0 -debug -p -g -traceback -check all -check bounds -fp-stack-check -assume realloc_lhs")
    set (COMPILE_FLAGS_HDF5toTEC_RELEASE "-O3 -assume realloc_lhs")
    set (LINK_FLAGS_HDF5toTEC_DEBUG "-g -traceback -check all -check bounds")
    set (LINK_FLAGS_HDF5toTEC_RELEASE "")

elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

    set (COMPILE_FLAGS_HDF5toTEC_DEBUG "-O0 -fcray-pointer -ffree-form -ffree-line-length-none -cpp")
    set (COMPILE_FLAGS_HDF5toTEC_RELEASE "-O3 -fcray-pointer -ffree-form -ffree-line-length-none -cpp")
    set (LINK_FLAGS_HDF5toTEC_DEBUG "")
    set (LINK_FLAGS_HDF5toTEC_RELEASE "")
endif()



# Need to re-include directory here
include_directories(${CMAKE_SOURCE_DIR}/src/tools)    # Need include here for macro-enhanced message handler messenger.h





add_executable(HDF5toTEC HDF5toTEC.f90)
target_link_libraries(HDF5toTEC coredg ${HDF5_LIBRARIES} ${TECPLOT_LIBRARIES})
set_target_properties(HDF5toTEC PROPERTIES 
                            COMPILE_FLAGS "${COMPILE_FLAGS_HDF5toTEC_RELEASE}"
                            LINK_FLAGS    "${LINK_FLAGS_HDF5toTEC_RELEASE} -lstdc++")





















