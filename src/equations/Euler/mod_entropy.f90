module mod_entropy
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: TWO, ZERO, THREE

    use type_chidg_data,    only: chidg_data_t

    use mod_interpolate,    only: interpolate_element
    use DNAD_D

    implicit none




contains



    !function compute_entropy_error(domain) result(error_sum)
    function compute_entropy_error(data) result(error_sum)
        type(chidg_data_t), intent(inout)   :: data


        !& DEBUG HARDCODED NUMBER OF NODES
        real(rk), dimension(data%mesh(1)%elems(1)%gq%vol%nnodes)  :: &
            rho, rhou, rhov, rhow, rhoE, p, entropy

        real(rk), dimension(data%mesh(1)%elems(1)%gq%vol%nnodes)    :: entropy_rise

        integer(ik) :: irho, irhou, irhov, irhow, irhoE
        integer(ik) :: ielem, iface, nelem, iseed, idom
        real(rk)    :: pinf, tinf, rhoinf, gam, entropy_ref, error_sum, vol, error, entropy_error, vol_sum


        pinf   = 110000._rk
        tinf   = 300._rk
        rhoinf = pinf/(tinf*287.15_rk)

        !& DEBUG - HARDCODED GAMMA
        gam    = 1.4_rk

        entropy_ref = pinf/(rhoinf**gam)


        !& DEBUG - HARDCODED EQUATION SET
        associate (mesh => data%mesh, sdata => data%sdata, eqnset => data%eqnset, prop => data%eqnset(1)%item%prop)


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
            do idom = 1,data%ndomains()

                nelem = data%mesh(idom)%nelem
                do ielem = 1,nelem

                    !
                    ! Interpolate variables to GQ nodes
                    !
                    call interpolate_element(mesh,sdata%q,idom,ielem,irho,  rho)
                    call interpolate_element(mesh,sdata%q,idom,ielem,irhou, rhou)
                    call interpolate_element(mesh,sdata%q,idom,ielem,irhov, rhov)
                    call interpolate_element(mesh,sdata%q,idom,ielem,irhow, rhow)
                    call interpolate_element(mesh,sdata%q,idom,ielem,irhoE, rhoE)


                    !
                    ! Compute pressure
                    !
                    call prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)


                    !
                    ! Compute entropy and entropy rise.
                    !
                    entropy = p/(rho**gam)
                    entropy_rise = ((entropy - entropy_ref)/entropy_ref)**TWO


                    !
                    ! Integrate entropy error
                    !
                    error = sum(entropy_rise * mesh(idom)%elems(ielem)%jinv * mesh(idom)%elems(ielem)%gq%vol%weights)


                    !
                    ! Compute element volume
                    !
                    vol = abs(sum(mesh(idom)%elems(ielem)%jinv * mesh(idom)%elems(ielem)%gq%vol%weights))


                    error_sum = error_sum + error
                    vol_sum   = vol_sum + vol


                end do ! ielem

            end do ! idom


        end associate


        error_sum = sqrt(error_sum/vol_sum)



    end function






end module mod_entropy
