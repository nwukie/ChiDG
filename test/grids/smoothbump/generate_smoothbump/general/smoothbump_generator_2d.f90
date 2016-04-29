program smoothbump_generator
    implicit none


    integer(kind=4)                 :: npt_x, npt_y, npt_z, &
                                       nelem_quartic_x, nelem_quartic_y, nelem_quartic_z, &
                                       ierr, ncoords, i, j, k, n, &
                                       ipt_x, ipt_y, ipt_z
    real(kind=8), allocatable       :: coord(:,:,:,:)
    real(kind=8)                    :: x,y,z,alpha
    character(len=100)              :: buffer



    !
    ! Get grid dimensions from user-input arguments
    !
    print*, 'Enter number of elements in each coordinate direction: nelem_x, nelem_y, nelem_z'
    read(*,*) nelem_quartic_x, nelem_quartic_y, nelem_quartic_z



    !
    ! Define number of quartic elements in each direction
    !
    !nelem_quartic_x = 6
    !nelem_quartic_y = 2
    !nelem_quartic_z = 2

    !
    ! Compute number of points in each direction
    !
    npt_x = nelem_quartic_x*4 + 1
    npt_y = nelem_quartic_y*4 + 1
    npt_z = nelem_quartic_z*4 + 1 


    !
    ! Allocate coordinate array
    !
    ncoords = 3
    allocate(coord(npt_x,npt_y,npt_z,ncoords),stat=ierr)
    if (ierr /= 0) stop "Allocation Error"


    !
    ! Create gridfile name
    !
    write(buffer,'( "smoothbump_",I0.0,"x",I0.0,"x",I0.0,".x")') nelem_quartic_x, nelem_quartic_y, nelem_quartic_z



    !
    ! Generate coordinates
    !
    do ipt_z = 1,npt_z
        
        do ipt_y = 1,npt_y




            do ipt_x = 1,npt_x

                x = -1.5_8 + real(ipt_x-1,kind=8)*(3._8 / real(npt_x-1,kind=8))
                alpha = 0.0625_8 * exp(-25._8 * (x**2._8))
                y = alpha + real(ipt_y-1,kind=8)*(0.8_8 - alpha)/real(npt_y-1,kind=8)
                z = 0._8 + real(ipt_z-1,kind=8)*(1._8 / real(npt_z-1,kind=8))


                coord(ipt_x,ipt_y,ipt_z,1) = x
                coord(ipt_x,ipt_y,ipt_z,2) = y
                coord(ipt_x,ipt_y,ipt_z,3) = z

            end do

        end do

    end do








    !
    ! Write coordinates to file
    !
    open(unit=7, file=trim(buffer), form='unformatted')

    write(7) 1
    write(7) npt_x, npt_y, npt_z
    write(7) (((( coord(i,j,k,n), i=1,npt_x), j=1,npt_y), k=1,npt_z), n=1,ncoords)

    close(7)






end program
