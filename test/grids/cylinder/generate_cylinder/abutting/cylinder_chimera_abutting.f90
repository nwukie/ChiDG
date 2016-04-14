program cylinder_chimera
    implicit none


    integer(kind=4), parameter  :: rk = 8
    real(rk),        parameter  :: PI = 3.14159265358979323846_rk

    integer(kind=4)                 :: npt_x, npt_y, npt_z, &
                                       nelem_quartic_xi, nelem_quartic_eta, nelem_quartic_zeta, &
                                       ierr, ncoords, i, j, k, n, &
                                       ipt_x, ipt_y, ipt_z
    real(rk), allocatable           :: coord(:,:,:,:), rotated_coord(:,:,:,:)
    real(rk)                        :: x,y,z,alpha
    character(len=100)              :: buffer

    real(rk)    :: xcenter, ycenter, zcenter, cradius, &
                   xmin, xmax, ymin, ymax, zmin, zmax, &
                   thetamin, thetamax, theta, clustering_factor

    real(rk), allocatable, dimension(:)   :: xupper, yupper, xlower, ylower



    xcenter = 0._rk
    ycenter = 0._rk
    zcenter = 0._rk

    cradius = 0.5_rk

    xmin = -3._rk
    xmax = 3._rk
    ymin = -3._rk
    ymax = 3._rk
    zmin = 0._rk
    zmax = 1._rk




    !
    ! Get grid dimensions from user-input arguments
    !
    print*, 'Enter number of elements in each coordinate direction: nelem_x, nelem_y, nelem_z'
    read(*,*) nelem_quartic_xi, nelem_quartic_eta, nelem_quartic_zeta



    !
    ! Compute number of points in each direction
    !
    npt_x = nelem_quartic_xi*4   + 1
    npt_y = nelem_quartic_eta*4  + 1
    npt_z = nelem_quartic_zeta*4 + 1 


    !
    ! Allocate coordinate array
    !
    ncoords = 3
    allocate(coord(npt_x,npt_y,npt_z,ncoords),stat=ierr)
    allocate(rotated_coord(npt_x,npt_y,npt_z,ncoords),stat=ierr)
    if (ierr /= 0) stop "Allocation Error"


    !
    ! Create gridfile name
    !
    write(buffer,'( "cylinder_chimera_abutting.x" )')



    !
    ! Generate upper and lower curves
    !
    thetamin = ((3._rk*PI/2._rk) + PI ) / 2._rk
    thetamax = ((3._rk*PI/2._rk) + 2*PI ) / 2._rk
    allocate(xupper(npt_x), yupper(npt_x), xlower(npt_x), ylower(npt_x))
    do ipt_x = 1,npt_x
        xlower(ipt_x) = xmin + real(ipt_x-1,rk)*((xmax-xmin) / real(npt_x-1,rk))
        ylower(ipt_x) = ymin

        theta = thetamin + real(ipt_x-1,rk)*( (thetamax - thetamin) / real(npt_x-1,rk))
        xupper(ipt_x) = cradius * dcos(theta)
        yupper(ipt_x) = cradius * dsin(theta)
    end do



    !
    ! Generate coordinates
    !
    do ipt_z = 1,npt_z
        
        do ipt_y = 1,npt_y




            do ipt_x = 1,npt_x

                clustering_factor = (real(ipt_y,rk)/real(npt_y,rk)) ** (-0.5_rk)

                x = xlower(ipt_x) + (clustering_factor)*real(ipt_y-1,rk)*((xupper(ipt_x)-xlower(ipt_x)) / real(npt_y-1,rk))
                y = ylower(ipt_x) + (clustering_factor)*real(ipt_y-1,rk)*((yupper(ipt_x)-ylower(ipt_x)) / real(npt_y-1,rk))
                z = 0._rk + real(ipt_z-1,rk)*(1._rk / real(npt_z-1,rk))


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

    !write(7) 1
    !write(7) npt_x, npt_y, npt_z
    write(7) 4
    write(7) npt_x, npt_y, npt_z, npt_x, npt_y, npt_z, npt_x, npt_y, npt_z, npt_x, npt_y, npt_z 

    !
    ! Write original points
    !
    write(7) (((( coord(i,j,k,n), i=1,npt_x), j=1,npt_y), k=1,npt_z), n=1,ncoords)

    !
    ! Perform 90-degree rotation
    !
    do ipt_z = 1,npt_z
        do ipt_y = 1,npt_y
            do ipt_x = 1,npt_x
                rotated_coord(ipt_x,ipt_y,ipt_z,1) = coord(ipt_x,ipt_y,ipt_z,1)*dcos(PI/2._rk) - &
                                                     coord(ipt_x,ipt_y,ipt_z,2)*dsin(PI/2._rk)

                rotated_coord(ipt_x,ipt_y,ipt_z,2) = coord(ipt_x,ipt_y,ipt_z,1)*dsin(PI/2._rk) + &
                                                     coord(ipt_x,ipt_y,ipt_z,2)*dcos(PI/2._rk)

                rotated_coord(ipt_x,ipt_y,ipt_z,3) = coord(ipt_x,ipt_y,ipt_z,3)
            end do
        end do
    end do
    write(7) (((( rotated_coord(i,j,k,n), i=1,npt_x), j=1,npt_y), k=1,npt_z), n=1,ncoords)





    !
    ! Perform 90-degree rotation
    !
    do ipt_z = 1,npt_z
        do ipt_y = 1,npt_y
            do ipt_x = 1,npt_x
                coord(ipt_x,ipt_y,ipt_z,1) = rotated_coord(ipt_x,ipt_y,ipt_z,1)*dcos(PI/2._rk) - &
                                             rotated_coord(ipt_x,ipt_y,ipt_z,2)*dsin(PI/2._rk)

                coord(ipt_x,ipt_y,ipt_z,2) = rotated_coord(ipt_x,ipt_y,ipt_z,1)*dsin(PI/2._rk) + &
                                             rotated_coord(ipt_x,ipt_y,ipt_z,2)*dcos(PI/2._rk)

                coord(ipt_x,ipt_y,ipt_z,3) = rotated_coord(ipt_x,ipt_y,ipt_z,3)
            end do
        end do
    end do
    write(7) (((( coord(i,j,k,n), i=1,npt_x), j=1,npt_y), k=1,npt_z), n=1,ncoords)



    !
    ! Perform 90-degree rotation
    !
    do ipt_z = 1,npt_z
        do ipt_y = 1,npt_y
            do ipt_x = 1,npt_x
                rotated_coord(ipt_x,ipt_y,ipt_z,1) = coord(ipt_x,ipt_y,ipt_z,1)*dcos(PI/2._rk) - &
                                                     coord(ipt_x,ipt_y,ipt_z,2)*dsin(PI/2._rk)

                rotated_coord(ipt_x,ipt_y,ipt_z,2) = coord(ipt_x,ipt_y,ipt_z,1)*dsin(PI/2._rk) + &
                                                     coord(ipt_x,ipt_y,ipt_z,2)*dcos(PI/2._rk)

                rotated_coord(ipt_x,ipt_y,ipt_z,3) = coord(ipt_x,ipt_y,ipt_z,3)
            end do
        end do
    end do
    write(7) (((( rotated_coord(i,j,k,n), i=1,npt_x), j=1,npt_y), k=1,npt_z), n=1,ncoords)
















    close(7)






end program cylinder_chimera
