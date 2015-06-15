module cgns
    use mod_kinds,  only: rk,ik
    implicit none

    include 'cgnslib_f.h'

contains


    subroutine plot3d_to_dgcgns()
      integer :: ni,nj,nk,k,j,i,index_file,ier,icelldim,iphysdim,index_base,isize,index_zone,index_coord
      real*8 x(21,17,9),y(21,17,9),z(21,17,9)
      dimension isize(3,3)
      character basename*32,zonename*32
!
!   create gridpoints for simple example:
      ni=21
      nj=17
      nk=9
      do k=1,nk
        do j=1,nj
          do i=1,ni
            x(i,j,k)=float(i-1)
            y(i,j,k)=float(j-1)
            z(i,j,k)=float(k-1)
          enddo
        enddo
      enddo
      write(6,'('' created simple 3-D grid points'')')
!
!   WRITE X, Y, Z GRID POINTS TO CGNS FILE
!   open CGNS file for write
      call cg_open_f('grid.cgns',CG_MODE_WRITE,index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
!   create base (user can give any name)
      basename='Base'
      icelldim=3
      iphysdim=3
      call cg_base_write_f(index_file,basename,icelldim,iphysdim,index_base,ier)
!   define zone name (user can give any name)
      zonename = 'Zone  1'
!   vertex size
      isize(1,1)=21
      isize(2,1)=17
      isize(3,1)=9
!   cell size
      isize(1,2)=isize(1,1)-1
      isize(2,2)=isize(2,1)-1
      isize(3,2)=isize(3,1)-1
!   boundary vertex size (always zero for structured grids)
      isize(1,3)=0
      isize(2,3)=0
      isize(3,3)=0
!   create zone
      call cg_zone_write_f(index_file,index_base,zonename,isize,Structured,index_zone,ier)
!   write grid coordinates (user must use SIDS-standard names here)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateX',x,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateY',y,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateZ',z,index_coord,ier)
!   close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote grid to file grid.cgns'')')
      stop




    end subroutine







end module cgns
