module mod_entropy
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, HALF, ONE, TWO, THREE, CYLINDRICAL
    use mod_interpolate,    only: interpolate_element_standard
    use mod_chidg_mpi,      only: ChiDG_COMM
    use mpi_f08,            only: MPI_AllReduce, MPI_REAL8, MPI_SUM
    use DNAD_D

    use type_chidg_data,    only: chidg_data_t
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


        real(rk), allocatable, dimension(:) :: &
            rho, rhou, rhov, rhow, rhoE, p, entropy

        !real(rk), dimension(data%mesh%domain(1)%elems(1)%gq%vol%nnodes)    :: entropy_rise
        real(rk), allocatable :: entropy_rise(:)

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
            irho  = prop%get_primary_field_index("Density"   )
            irhou = prop%get_primary_field_index("Momentum-1")
            irhov = prop%get_primary_field_index("Momentum-2")
            irhow = prop%get_primary_field_index("Momentum-3")
            irhoE = prop%get_primary_field_index("Energy"    )


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
            do idom = 1,data%mesh%ndomains()

                nelem = data%mesh%domain(idom)%nelem
                do ielem = 1,nelem

                    !
                    ! Interpolate variables to GQ nodes
                    !
                    rho  = interpolate_element_standard(mesh,sdata%q,idom,ielem,irho, 1, 'value')
                    rhou = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhou,1, 'value')
                    rhov = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhov,1, 'value')
                    rhow = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhow,1, 'value')
                    rhoE = interpolate_element_standard(mesh,sdata%q,idom,ielem,irhoE,1, 'value')

                    if (mesh%domain(idom)%elems(ielem)%coordinate_system == CYLINDRICAL) then
                        rhov = rhov / mesh%domain(idom)%elems(ielem)%interp_coords_def(:,1)
                    end if


                    !
                    ! Compute pressure
                    !
                    !p = prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)
                    p = (1.4_rk-ONE)*(rhoE - HALF*( (rhou*rhou) + (rhov*rhov) + (rhow*rhow) )/rho )


                    !
                    ! Compute entropy and entropy rise.
                    !
                    entropy = p/(rho**gam)
                    entropy_rise = ((entropy - entropy_ref)/entropy_ref)**TWO


                    !
                    ! Integrate entropy error
                    !
                    error = sum(entropy_rise * mesh%domain(idom)%elems(ielem)%jinv_def * mesh%domain(idom)%elems(ielem)%basis_s%weights_element())


                    !
                    ! Compute element volume
                    !
                    vol = abs(sum(mesh%domain(idom)%elems(ielem)%jinv_def * mesh%domain(idom)%elems(ielem)%basis_s%weights_element()))


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
