module type_fluid_pseudo_timestep
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: THIRD
    use type_pseudo_timestep,   only: pseudo_timestep_t
    use type_mesh,              only: mesh_t
    use type_properties,        only: properties_t
    use type_solverdata,        only: solverdata_t
    use mod_interpolate,        only: interpolate_element_standard
    implicit none



    !>  Pseudo Time-step calculator for fluid equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    type, extends(pseudo_timestep_t), public :: fluid_pseudo_timestep_t

    contains

        procedure   :: compute

    end type fluid_pseudo_timestep_t
    !***********************************************************************************





contains


    !> Routine to compute the local time-step in each element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[inout]   domain      domain_t instance containing mesh and solution data
    !!
    !-----------------------------------------------------------------------------------
    subroutine compute(self,idomain,mesh,prop,sdata,cfl,itime)
        class(fluid_pseudo_timestep_t), intent(in)      :: self
        integer(ik),                    intent(in)      :: idomain
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(properties_t),             intent(in)      :: prop
        type(solverdata_t),             intent(inout)   :: sdata
        real(rk),                       intent(in)      :: cfl
        integer(ik),                    intent(in)      :: itime


        integer(ik) :: ielem

        integer(ik) :: irho, irhou, irhov, irhow, irhoE

        real(rk), allocatable, dimension(:) ::  &
                rho, rhou, rhov, rhow, rhoE,    &
                gam, p, vmag, c


        real(rk)    ::  h, &    !< element spacing parameter
                        lam     !< characteristic speed

        
        !
        ! Get variable indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("X-Momentum")
        irhov = prop%get_primary_field_index("Y-Momentum")
        irhow = prop%get_primary_field_index("Z-Momentum")
        irhoE = prop%get_primary_field_index("Energy"    )



       do ielem = 1,mesh(idomain)%nelem


            !
            ! Interpolate variables
            !
            rho  = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irho, itime, 'value')
            rhou = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhou,itime, 'value')
            rhov = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhov,itime, 'value')
            rhow = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhow,itime, 'value')
            rhoE = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhoE,itime, 'value')


            !
            ! Compute pressure
            !
            p = prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)
        

            !
            ! Compute cell sound speed
            !
            gam = prop%fluid%compute_gamma(rho,rhou,rhov,rhow,rhoE)


            

            ! Compiling with DEBUG and bounds checking, gfortran will say 'c' is not correct size.
            ! This is not correct because 'c' should be sized according to the rhs of the expression.
            ! The sizes of gam, p, and rho are all the same. This is a recognized bug.
            !
            !   GCC/GFortran Bugzilla Bug 52162 
            !
            !   It gets triggered by calling the intrinsic sqrt. 
            !
            ! The explicit 'allocate's and 'deallocate's here work around the bug, but really
            ! should not be needed, since the arrays are automatically allocated. The compiled
            ! code just doesn't recognize the size of the automatically allocated arrays
            ! correctly.
            !
            if (allocated(c)) deallocate(c)
            if (allocated(vmag)) deallocate(vmag)
            allocate(c(size(p)), vmag(size(p)))
            c = sqrt(gam * p / rho)


            !
            ! Compute velocity magnitude
            !
            vmag = sqrt((rhou*rhou + rhov*rhov + rhow*rhow)/(rho*rho))


            !
            ! Compute mean characteristic speed. First compute average velocity magnitude and sound speed
            !
            lam = sum(vmag)/size(vmag) + sum(c)/size(vmag)


            !
            ! Compute element spacing parameter
            !
            h = mesh(idomain)%elems(ielem)%vol**(THIRD)


            !
            ! Compute elemen-local timestep
            !
            sdata%dt(idomain,ielem) = (cfl*h)/lam



        end do  ! ielem



    end subroutine compute
    !*************************************************************************************************













end module type_fluid_pseudo_timestep
