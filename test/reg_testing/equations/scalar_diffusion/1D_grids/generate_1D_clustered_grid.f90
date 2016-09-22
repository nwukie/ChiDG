program generate_1D_clustered_grid
    implicit none

    real(kind=8), parameter         :: PI  = 3.14159265358979323846_8
    real(kind=8), parameter         :: ONE = 1._8
    real(kind=8), parameter         :: TWO = 2._8


    integer(kind=4)                 :: npt_x, npt_y, npt_z, &
                                       nelem_x, nelem_y, nelem_z, &
                                       ierr, ncoords, i, j, k, n, &
                                       ipt_x, ipt_y, ipt_z, mapping
    real(kind=8), allocatable       :: coord(:,:,:,:)
    real(kind=8)                    :: x,y,z,alpha
    character(len=100)              :: buffer



    !
    ! Get grid dimensions from user-input arguments
    !
    print*, 'Enter number of elements in each coordinate direction: nelem_x, nelem_y, nelem_z'
    read(*,*) nelem_x, nelem_y, nelem_z


    !
    ! Compute number of points in each direction
    !
    mapping = 1 ! linear

    npt_x = nelem_x*mapping + 1
    npt_y = nelem_y*mapping + 1
    npt_z = nelem_z*mapping + 1 


    !
    ! Allocate coordinate array
    !
    ncoords = 3
    allocate(coord(npt_x,npt_y,npt_z,ncoords),stat=ierr)
    if (ierr /= 0) stop "Allocation Error"


    !
    ! Create gridfile name
    !
    write(buffer,'( "1D_clustered_grid.x")')



    !
    ! Generate coordinates
    !
    do ipt_z = 1,npt_z
        do ipt_y = 1,npt_y
            do ipt_x = 1,npt_x

!                x = 0._8 + real(ipt_x-1,kind=8)*(1._8   / real(npt_x-1,kind=8))
!                if (ipt_x == 1) then
!                    x = 0._8
!                else if (ipt_x == npt_x) then
!                    x = 1._8
!                else
                    !x = 1._8 - tanh((PI/TWO)*(ONE - (real(ipt_x,kind=8)/(real(npt_x,kind=8)-ONE))))/tanh(PI/TWO)
                    x = 1._8 - tanh((PI/TWO)*(ONE - (real(ipt_x-1,kind=8)/real(npt_x-1,kind=8)) ))/tanh(PI/TWO)
!                end if
                y = 0._8 + real(ipt_y-1,kind=8)*(0.1_8 / real(npt_y-1,kind=8))
                z = 0._8 + real(ipt_z-1,kind=8)*(0.1_8 / real(npt_z-1,kind=8))


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






end program generate_1D_clustered_grid
