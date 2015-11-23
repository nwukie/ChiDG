module mod_partitioners
    use type_meshdata,  only: meshdata_t
    use mod_mpichi





contains




    !>
    !!
    !!
    !!
    !! 
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine partition_mesh(meshdata)
        type(meshdata_t),   intent(inout)   :: meshdata(:)

        integer :: ndomains

        logical :: serial = .false.
        logical :: parallel = .false.

        logical :: partition_blocks = .false.

        type(meshdata_t), allocatable   :: meshdata_tmpA(:)
        type(meshdata_t), allocatable   :: meshdata_tmpB(:)
        type(meshdata_t), allocatable   :: meshdata_split(:)




        !
        ! Block splitting with master process
        !
        if ( irank == GLOBAL_MASTER ) then


            !
            ! Get original number of block domains
            !
            ndomains = size(meshdata)


            !
            ! Is this a serial or parallel calculation
            !
            if ( nrank > 1 ) then
                parallel = .true.
                serial   = .false.
            else
                serial   = .true.
                parallel = .false.
            end if



            !
            ! Check if we need to split blocks
            !
            if ( ndomains > nrank ) partition_blocks = .true.




            !
            ! Action to split blocks
            !
            if (partition_blocks) then



            else



            end if




        end if ! irank == GLOBAL_MASTER




















        !
        ! Block send/recv
        !











    end subroutine
    !----------------------------------------------------------------------------------------------














    subroutine assign_processors(meshdata)
        type(meshdata_t),   intent(in)  :: meshdata(:)

        integer :: ndomains, idom

        ndomains = size(meshdata)


        if ( ndomains == nrank ) then

            !
            ! Assign consecutive processor numbers
            !
            do idom = 1,ndomains
                meshdata(idom)%proc = idom
            end do

        else

            
            !
            ! Agglomerate blocks on some processors
            !




            

        end if


    end subroutine assign_processors



























end module mod_partitioners
