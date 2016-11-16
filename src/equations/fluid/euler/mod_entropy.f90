module mod_entropy
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: TWO, ZERO, THREE

    use type_chidg_data,    only: chidg_data_t

    use mod_interpolate,    only: interpolate_element_standard
    use mod_chidg_mpi,      only: ChiDG_COMM
    use mpi_f08,            only: MPI_AllReduce, MPI_REAL8, MPI_SUM
    use DNAD_D
    implicit none




contains



    !> Function to compute the entropy error in a domain. Assuming the Euler equations
    !! and isentropic flow.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    function compute_entropy_error(data) result(entropy_error)
        type(chidg_data_t), intent(inout)   :: data


        !& DEBUG HARDCODED NUMBER OF NODES
        !real(rk), dimension(data%mesh(1)%elems(1)%gq%vol%nnodes)  :: &
        real(rk), allocatable, dimension(:) :: &
            rho, rhou, rhov, rhow, rhoE, p, entropy

        real(rk), dimension(data%mesh(1)%elems(1)%gq%vol%nnodes)    :: entropy_rise

        integer(ik) :: irho, irhou, irhov, irhow, irhoE
        integer(ik) :: ielem, iface, nelem, iseed, idom, ierr
        real(rk)    :: pinf, tinf, rhoinf, gam, entropy_ref, error_sum, error_sum_reduced, vol, error, entropy_error, vol_sum, vol_sum_reduced


        pinf   = 110000._rk
        tinf   = 300._rk
        rhoinf = pinf/(tinf*287.15_rk)

        !& DEBUG - HARDCODED GAMMA
        gam    = 1.4_rk

        entropy_ref = pinf/(rhoinf**gam)


        !& DEBUG - HARDCODED EQUATION SET
        associate (mesh => data%mesh, sdata => data%sdata, eqnset => data%eqnset, prop => data%eqnset(1)%prop)


            ! Get equation indices
            irho  = prop%get_equation_index("Density"   )
            irhou = prop%get_equation_index("X-Momentum")
            irhov = prop%get_equation_index("Y-Momentum")
            irhow = prop%get_equation_index("Z-Momentum")
            irhoE = prop%get_equation_index("Energy"    )


            !
            ! No need to seed element derivatives
            !
            iseed = 0




            !
            ! Zero entropy error
            !
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
                    rho  = interpolate_element_standard(mesh,sdata%q,idom,ielem,irho, 1,  'value')
                    rhou = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhou,1, 'value')
                    rhov = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhov,1, 'value')
                    rhow = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhow,1, 'value')
                    rhoE = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhoE,1, 'value')


                    !
                    ! Compute pressure
                    !
                    p = prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)


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



        ! Reduce the total error across processors
        call MPI_AllReduce(error_sum,error_sum_reduced,1,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)

        ! Reduce the total volume across processors
        call MPI_AllReduce(vol_sum,vol_sum_reduced,1,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)

        ! Compute the global volume-weighted entropy error
        entropy_error = sqrt(error_sum_reduced/vol_sum_reduced)


    end function compute_entropy_error
    !**********************************************************************************************






end module mod_entropy
