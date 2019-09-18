module mod_modify_jacobian_sst
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use mod_hdfio,              only: write_fields_hdf
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER, IRANK, NRANK
    use mod_io,                 only: verbosity
    use mpi_f08,                only: MPI_Barrier

    use type_chidg_data,        only: chidg_data_t
    use type_chidg_vector

    use ieee_arithmetic,        only: ieee_is_nan
    implicit none

contains

    subroutine modify_jacobian_sst(data)
        class(chidg_data_t),                     intent(inout)           :: data

        integer(ik) :: idom, eqn_ID, ielem, itime, imat, ieqn, ieqn_col, ieqn_sst, nterms, rstart, rend, cstart, cend
        character(len=24), dimension(2) :: sst_fields = [character(len=24) :: 'Density * Omega',  &
                        'Density * k']

        character(100) :: field_name
        logical :: row_is_sst, row_is_meanflow, col_is_sst, col_is_meanflow, decouple

        ! Decouple Mean Flow and sst Jacobians 
        do idom = 1,data%mesh%ndomains()
            do ielem = 1,data%mesh%domain(idom)%nelem
                eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                do itime = 1,data%mesh%domain(idom)%ntime
                    do imat = 1,data%sdata%lhs%dom(idom)%lblks(ielem,itime)%size()
                        do ieqn = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()
                            field_name = data%eqnset(eqn_ID)%prop%get_primary_field_name(ieqn)
                            row_is_sst = .false.
                            do ieqn_sst = 1,size(sst_fields)
                               if (trim(field_name) == trim(sst_fields(ieqn_sst))) row_is_sst = .true.
                            end do

                            row_is_meanflow = (.not. row_is_sst)
                            ! Need to compute row and column extends in diagonal so we can
                            ! selectively apply the mass matrix to the sub-block diagonal
                            nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s
                            rstart = 1 + (ieqn-1) * nterms
                            rend   = (rstart-1) + nterms

                            do ieqn_col = 1, data%eqnset(eqn_ID)%prop%nprimary_fields()
                                field_name = data%eqnset(eqn_ID)%prop%get_primary_field_name(ieqn_col)
                                col_is_sst = .false.
                                do ieqn_sst = 1,size(sst_fields)
                                   if (trim(field_name) == trim(sst_fields(ieqn_sst))) col_is_sst = .true.
                                end do

                                col_is_meanflow = (.not. col_is_sst)

                                cstart = 1 + (ieqn_col-1) * nterms
                                cend   = (cstart-1) + nterms

                                ! Zero out MF-sst matrix entries
                                decouple = ((row_is_meanflow .and. col_is_sst) .or. (row_is_sst .and. col_is_meanflow))
                                !decouple = (row_is_meanflow .and. col_is_sst)
                                if (decouple) then
                                    data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  = ZERO
                                end if

                            end do

                        end do !ieqn
                    end do !imat
                end do !itime
            end do !ielem
        end do !idom


    end subroutine
end module mod_modify_jacobian_sst
