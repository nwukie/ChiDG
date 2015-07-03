add_lib_sources(hdf5_dgtools.f90 mod_io.f90)



add_executable(plot3d_to_hdf5 plot3d_to_hdf5.f90 ../constants/mod_kinds.f90)
target_link_libraries(plot3d_to_hdf5 ${HDF5_LIBRARIES})
set_target_properties(plot3d_to_hdf5 PROPERTIES 
                            COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}"
                            LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")



