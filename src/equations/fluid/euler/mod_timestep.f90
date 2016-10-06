module mod_timestep
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: THIRD
    use type_chidg_data,    only: chidg_data_t
    use mod_interpolate,    only: interpolate_element_standard
    implicit none



contains


    !> Routine to compute the local time-step in each element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[inout]   domain      domain_t instance containing mesh and solution data
    !!
    !-----------------------------------------------------------------------------------
    subroutine compute_timestep(data,cfl)
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           intent(in)      :: cfl


        integer(ik)     :: ielem, nelem, idom

        integer(ik) :: irho, irhou, irhov, irhow, irhoE

        real(rk), allocatable, dimension(:) ::  &
                rho, rhou, rhov, rhow, rhoE,    &
                gam, p, vmag, c


        real(rk)    ::  h, &    !< element spacing parameter
                        lam     !< characteristic speed

        
        !
        ! Loop through elements and compute time-step function
        !
        do idom = 1,data%ndomains()

            !
            ! Get variable indices
            !
            irho  = data%eqnset(idom)%prop%get_equation_index("Density"   )
            irhou = data%eqnset(idom)%prop%get_equation_index("X-Momentum")
            irhov = data%eqnset(idom)%prop%get_equation_index("Y-Momentum")
            irhow = data%eqnset(idom)%prop%get_equation_index("Z-Momentum")
            irhoE = data%eqnset(idom)%prop%get_equation_index("Energy"    )



           nelem = data%mesh(idom)%nelem
           do ielem = 1,nelem


                !
                ! Interpolate variables
                !
                rho  = interpolate_element_standard(data%mesh,data%sdata%q,idom,ielem,irho,  'value')
                rhou = interpolate_element_standard(data%mesh,data%sdata%q,idom,ielem,irhou, 'value')
                rhov = interpolate_element_standard(data%mesh,data%sdata%q,idom,ielem,irhov, 'value')
                rhow = interpolate_element_standard(data%mesh,data%sdata%q,idom,ielem,irhow, 'value')
                rhoE = interpolate_element_standard(data%mesh,data%sdata%q,idom,ielem,irhoE, 'value')


                !
                ! Compute pressure
                !
                p = data%eqnset(idom)%prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)
            

                !
                ! Compute cell sound speed
                !
                gam = data%eqnset(idom)%prop%fluid%compute_gamma(rho,rhou,rhov,rhow,rhoE)


                

                ! Compiling with DEBUG and bounds checking, gfortran will say 'c' is not correct size.
                ! This is not correct because 'c' should be sized according to the rhs of the expression.
                ! The sizes of gam, p, and rho are all the same. This is a recognized bug.
                !
                !   GCC/GFortran Bugzilla Bug 52162 
                !
                !   It gets triggered by calling the intrinsic sqrt. 
                !
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
                h = data%mesh(idom)%elems(ielem)%vol**(THIRD)


                !
                ! Compute elemen-local timestep
                !
                data%sdata%dt(idom,ielem) = (cfl*h)/lam
                !data%sdata%dt(idom,ielem) = cfl*h



            end do  ! ielem
        end do  ! idom



    end subroutine compute_timestep
    !*************************************************************************************************













end module mod_timestep
