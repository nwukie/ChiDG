module mod_timestep
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: THIRD
    use type_domain,        only: domain_t

    use mod_interpolate,    only: interpolate

    implicit none



contains


    !> Routine to compute the local time-step in each element
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   domain      domain_t instance containing mesh and solution data
    !!
    !-----------------------------------------------------------------------------------
    subroutine compute_timestep(domain,cfl)
        type(domain_t),     intent(inout)   :: domain
        real(rk),           intent(in)      :: cfl


        integer(ik) :: ielem, nelem

        integer(ik) :: irho, irhou, irhov, irhow, irhoE

        !& DEBUG: HARD CODED GQ NODES BASED ON FIRST ELEMENT
        real(rk), dimension(domain%mesh%elems(1)%gq%vol%nnodes)  :: rho, rhou, rhov, rhow, rhoE, &
                                                                    c,   &   !< mean sound speed
                                                                    gam, &   !< ratio of specific heats
                                                                    p,   &   !< pressure
                                                                    vmag     !< velocity magnitude

        real(rk)    ::  h, &    !< element spacing parameter
                        lam     !< characteristic speed

        !
        ! Get number of elements
        !
        nelem = domain%mesh%nelem
        

        !
        ! Get variable indices
        !
        irho  = domain%eqnset%prop%get_eqn_index("rho")
        irhou = domain%eqnset%prop%get_eqn_index("rhou")
        irhov = domain%eqnset%prop%get_eqn_index("rhov")
        irhow = domain%eqnset%prop%get_eqn_index("rhow")
        irhoE = domain%eqnset%prop%get_eqn_index("rhoE")


        associate ( elems => domain%mesh%elems, q => domain%sdata%q )
        !
        ! Loop through elements and compute time-step function
        !
        do ielem = 1,nelem


            !
            ! Interpolate variables
            !
            call interpolate(elems,q,ielem,irho,  rho)
            call interpolate(elems,q,ielem,irhou, rhou)
            call interpolate(elems,q,ielem,irhov, rhov)
            call interpolate(elems,q,ielem,irhow, rhow)
            call interpolate(elems,q,ielem,irhoE, rhoE)


            !
            ! Compute pressure
            !
            call domain%eqnset%prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)
        

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
            h = elems(ielem)%vol**(THIRD)


            !
            ! Compute elemen-local timestep
            !
            domain%sdata%dt(ielem) = (cfl*h)/lam



        end do
        end associate


    end subroutine













end module mod_timestep
