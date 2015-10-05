module mod_entropy
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: TWO, ZERO, THREE
    use type_domain,        only: domain_t

    use mod_interpolate,    only: interpolate
    use DNAD_D

    implicit none




contains



    function compute_entropy_error(domain) result(error_sum)
        type(domain_t), intent(inout)   :: domain


        !& DEBUG HARDCODED NUMBER OF NODES
        type(AD_D), dimension(domain%mesh%elems(1)%gq%vol%nnodes)  :: &
            rho, rhou, rhov, rhow, rhoE, p, entropy

        real(rk), dimension(domain%mesh%elems(1)%gq%vol%nnodes)    :: entropy_rise

        integer(ik) :: irho, irhou, irhov, irhow, irhoE
        integer(ik) :: ielem, iface, nelem, iseed
        real(rk)    :: pinf, tinf, rhoinf, gam, entropy_ref, error_sum, vol, error, entropy_error, vol_sum


        pinf   = 110000._rk
        tinf   = 300._rk
        rhoinf = pinf/(tinf*287.15_rk)

        !& DEBUG - HARDCODED GAMMA
        gam    = 1.4_rk

        entropy_ref = pinf/(rhoinf**gam)


        associate (mesh => domain%mesh, sdata => domain%sdata, eqnset => domain%eqnset, prop => domain%eqnset%prop)


            ! Get equation indices
            irho  = prop%get_eqn_index('rho')
            irhou = prop%get_eqn_index('rhou')
            irhov = prop%get_eqn_index('rhov')
            irhow = prop%get_eqn_index('rhow')
            irhoE = prop%get_eqn_index('rhoE')


            !
            ! No need to seed element derivatives
            !
            iseed = 0




            !
            ! Zero entropy error
            !
            entropy_error = ZERO
            error_sum     = ZERO
            vol           = ZERO
            vol_sum       = ZERO




            !
            ! Loop over elements and accumulate entropy error
            !
            nelem = mesh%nelem
            do ielem = 1,nelem

                !
                ! Interpolate variables to GQ nodes
                !
                call interpolate(mesh%elems,sdata%q,ielem,irho,  rho,  iseed)
                call interpolate(mesh%elems,sdata%q,ielem,irhou, rhou, iseed)
                call interpolate(mesh%elems,sdata%q,ielem,irhov, rhov, iseed)
                call interpolate(mesh%elems,sdata%q,ielem,irhow, rhow, iseed)
                call interpolate(mesh%elems,sdata%q,ielem,irhoE, rhoE, iseed)


                !
                ! Compute pressure
                !
                call prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)



                !
                ! Compute entropy and entropy rise.
                !
                entropy = p/(rho**gam)
                entropy_rise = ((entropy(:)%x_ad_ - entropy_ref)/entropy_ref)**TWO




                !
                ! Integrate entropy error
                !
                error = sum(entropy_rise * mesh%elems(ielem)%jinv * mesh%elems(ielem)%gq%vol%weights)



                !
                ! Compute element volume
                !
                vol = abs(sum(mesh%elems(ielem)%jinv * mesh%elems(ielem)%gq%vol%weights))


!                print*, 'Entropy', error
!                print*, 'Volume', vol

                error_sum = error_sum + error
                vol_sum   = vol_sum + vol



            end do


        end associate


        error_sum = sqrt(error_sum/vol_sum)



    end function






end module mod_entropy
