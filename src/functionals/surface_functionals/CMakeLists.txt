add_subdirectory(mass_averaged_functionals)
add_subdirectory(scalar_functionals)

add_lib_sources(
                fcl_mass_flux.f90
                fcl_mass_flux_balance.f90
                fcl_xforce.f90
                fcl_yforce.f90
                fcl_zforce.f90
                fcl_xforce_coeff.f90
                fcl_yforce_coeff.f90
                fcl_zforce_coeff.f90
                )
