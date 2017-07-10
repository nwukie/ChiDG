module type_fluid_pseudo_timestep
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: THIRD, ONE, HALF
    use mod_fluid,              only: omega, gam
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
        type(mesh_t),               intent(inout)   :: mesh
        type(properties_t),             intent(in)      :: prop
        type(solverdata_t),             intent(inout)   :: sdata
        real(rk),                       intent(in)      :: cfl(:)
        integer(ik),                    intent(in)      :: itime


        integer(ik) :: ielem, ieqn

        integer(ik) :: irho, irhou, irhov, irhow, irhoE

        real(rk), allocatable, dimension(:) ::  &
                rho, rhou, rhov, rhow, rhoE,    &
                p, vmag, c, r


        real(rk)    ::  h, lam

        
        !
        ! Get variable indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("Momentum-1")
        irhov = prop%get_primary_field_index("Momentum-2")
        irhow = prop%get_primary_field_index("Momentum-3")
        irhoE = prop%get_primary_field_index("Energy"    )



       do ielem = 1,mesh%domain(idomain)%nelem


            !
            ! Interpolate variables
            !
            rho  = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irho, itime, 'value')
            rhou = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhou,itime, 'value')
            rhov = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhov,itime, 'value')
            rhow = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhow,itime, 'value')
            rhoE = interpolate_element_standard(mesh,sdata%q,idomain,ielem,irhoE,itime, 'value')


            !
            ! Account for Cylindrical coordinates. Get tangential momentum from angular momentum
            !
            if ( mesh%domain(idomain)%elems(ielem)%coordinate_system == 'Cylindrical' ) then
                !rhov = rhov / mesh%domain(idomain)%elems(ielem)%quad_pts(:,1)%c1_
                rhov = rhov / mesh%domain(idomain)%elems(ielem)%quad_pts(:,1)
            end if

            !
            ! Compute pressure
            !
            p = (gam - ONE)*(rhoE - HALF*(rhou*rhou + rhov*rhov + rhow*rhow)/rho )

            

            ! Compiling with DEBUG and bounds checking, gfortran will say 'c' is not 
            ! correct size. This is not correct because 'c' should be sized according to 
            ! the rhs of the expression. The sizes of gam, p, and rho are all the same. 
            ! This is a recognized bug.
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
            !r = mesh%domain(idomain)%elems(ielem)%quad_pts(:)%c1_
            r = mesh%domain(idomain)%elems(ielem)%quad_pts(:,1)
            vmag = sqrt((rhou*rhou + rhov*rhov + rhow*rhow)/(rho*rho))
            !vmag = sqrt((rhou*rhou + (rhov-rho*omega*r)*(rhov-rho*omega*r) + rhow*rhow)/(rho*rho))


            !
            ! Compute mean characteristic speed. First compute average velocity 
            ! magnitude and sound speed
            !
            lam = sum(vmag)/size(vmag) + sum(c)/size(vmag)


            !
            ! Compute element spacing parameter
            !
            h = mesh%domain(idomain)%elems(ielem)%vol**(THIRD)


            !
            ! Compute elemen-local timestep
            !
            !sdata%dt(idomain,ielem) = (cfl*h)/lam
            do ieqn = 1,size(cfl)
                mesh%domain(idomain)%elems(ielem)%dtau(ieqn) = cfl(ieqn)*h/lam
            end do



        end do  ! ielem



    end subroutine compute
    !***************************************************************************************













end module type_fluid_pseudo_timestep
