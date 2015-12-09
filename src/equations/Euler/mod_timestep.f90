module mod_timestep
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: THIRD
    use type_chidg_data,    only: chidg_data_t
    use type_seed,          only: seed_t

    use mod_interpolate,    only: interpolate_element

    implicit none



contains


    !> Routine to compute the local time-step in each element
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   domain      domain_t instance containing mesh and solution data
    !!
    !-----------------------------------------------------------------------------------
    subroutine compute_timestep(data,cfl)
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           intent(in)      :: cfl


        integer(ik)     :: ielem, nelem, idom
        type(seed_t)    :: seed

        integer(ik) :: irho, irhou, irhov, irhow, irhoE

        !& DEBUG: HARD CODED GQ NODES BASED ON FIRST ELEMENT
        real(rk), dimension(data%mesh(1)%elems(1)%gq%vol%nnodes)  :: rho, rhou, rhov, rhow, rhoE, &
                                                                    c,   &   !< mean sound speed
                                                                    gam, &   !< ratio of specific heats
                                                                    p,   &   !< pressure
                                                                    vmag     !< velocity magnitude

        real(rk)    ::  h, &    !< element spacing parameter
                        lam     !< characteristic speed

        
        !associate ( elems => domain%mesh%elems, q => domain%sdata%q )
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
                seed%idom  = 0
                seed%ielem = 0
                call interpolate_element(data%mesh,data%sdata%q,idom,ielem,irho,  rho)
                call interpolate_element(data%mesh,data%sdata%q,idom,ielem,irhou, rhou)
                call interpolate_element(data%mesh,data%sdata%q,idom,ielem,irhov, rhov)
                call interpolate_element(data%mesh,data%sdata%q,idom,ielem,irhow, rhow)
                call interpolate_element(data%mesh,data%sdata%q,idom,ielem,irhoE, rhoE)


                !
                ! Compute pressure
                !
                call data%eqnset(idom)%item%prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)
            

                !
                ! Compute cell sound speed
                !
                !& DEBUG - HARDCODED GAMMA
                gam = 1.4_rk
                c = sqrt(gam * p / rho)


                !
                ! Compute velocity magnitude
                !
                vmag = sqrt((rhou*rhou + rhov*rhov + rhow*rhow)/(rho*rho))


                !
                ! Compute mean characteristic speed. First compute average velocity magnitude and sound speed
                !
                ! lam = vmag + c
                lam = sum(vmag)/size(vmag) + sum(c)/size(vmag)


                !
                ! Compute element spacing parameter
                !
                h = data%mesh(idom)%elems(ielem)%vol**(THIRD)


                !
                ! Compute elemen-local timestep
                !
                !domain%sdata%dt(ielem) = (cfl*h)/lam
                data%sdata%dt(idom,ielem) = (cfl*h)/lam



            end do  ! ielem
        end do  ! idom

        !end associate


    end subroutine













end module mod_timestep
