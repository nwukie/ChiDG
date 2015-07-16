add_lib_sources(mod_io.f90
                mod_hdfio.f90
                mod_tecio_interface.f90
                mod_tecio.f90)












#------------------------------------------------------------------                                                                        
#                     Set compiler flags for 'plot3d_to_hdf5'
#------------------------------------------------------------------                                                                        
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (COMPILE_FLAGS_plot3dtohdf5_DEBUG   "-O0 -debug -p -g -traceback -check all -check bounds -fp-stack-check -assume realloc_lhs")
    set (COMPILE_FLAGS_plot3dtohdf5_RELEASE "-O3 -assume realloc_lhs")
    set (LINK_FLAGS_plot3dtohdf5_DEBUG "-g -traceback -check all -check bounds")
    set (LINK_FLAGS_plot3dtohdf5_RELEASE "")

elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

    set (COMPILE_FLAGS_plot3dtohdf5_DEBUG "-O3 -fcray-pointer -ffree-form -ffree-line-length-none -cpp")
    set (COMPILE_FLAGS_plot3dtohdf5_RELEASE "")
    set (LINK_FLAGS_plot3dtohdf5_DEBUG "")
    set (LINK_FLAGS_plot3dtohdf5_RELEASE "")
endif()



add_executable(plot3d_to_hdf5 plot3d_to_hdf5.f90 ../constants/mod_kinds.f90)
target_link_libraries(plot3d_to_hdf5 ${HDF5_LIBRARIES})
set_target_properties(plot3d_to_hdf5 PROPERTIES 
                            COMPILE_FLAGS "${COMPILE_FLAGS_plot3dtohdf5_RELEASE}"
                            LINK_FLAGS    "${LINK_FLAGS_plot3dtohdf5_RELEASE}")



