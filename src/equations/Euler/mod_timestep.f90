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

        !& DEBUG: HARD CODED GQ NODES BASED ON FIRST ELEMENT
        !real(rk), dimension(data%mesh(1)%elems(1)%gq%vol%nnodes)  :: &
        real(rk), allocatable, dimension(:) ::  &
                rho, rhou, rhov, rhow, rhoE,    &
                c, gam, p, vmag, tmp

        real(rk)    ::  h, &    !< element spacing parameter
                        lam     !< characteristic speed

        
        !
        ! Loop through elements and compute time-step function
        !
        do idom = 1,data%ndomains()

        !
        ! Get variable indices
        !
        irho  = data%eqnset(idom)%item%prop%get_eqn_index("rho")
        irhou = data%eqnset(idom)%item%prop%get_eqn_index("rhou")
        irhov = data%eqnset(idom)%item%prop%get_eqn_index("rhov")
        irhow = data%eqnset(idom)%item%prop%get_eqn_index("rhow")
        irhoE = data%eqnset(idom)%item%prop%get_eqn_index("rhoE")



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
                call data%eqnset(idom)%item%prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)
            

                !
                ! Compute cell sound speed
                !
                call data%eqnset(idom)%item%prop%fluid%compute_gamma(rho,rhou,rhov,rhow,rhoE,gam)



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



            end do  ! ielem
        end do  ! idom



    end subroutine compute_timestep
    !*************************************************************************************************













end module mod_timestep
