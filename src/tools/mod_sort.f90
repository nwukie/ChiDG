module mod_sort
    use mod_types,  only: rk,ik
    implicit none

contains
! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      real(kind=rk), DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      real(kind=rk)                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)       ! assume the first is the min
      Location = Start          ! record its position
      DO i = Start+1, End       ! start with next elements
         IF (x(i) < Minimum) THEN   !   if x(i) less than the min?
            Minimum  = x(i)     !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real(kind=rk), INTENT(INOUT) :: a, b
      real(kind=rk)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x,y,Size)
      IMPLICIT  NONE
      real(kind=rk), DIMENSION(1:), INTENT(INOUT) :: x,y
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1          ! except for the last
         Location = FindMinimum(x, i, Size) ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
         CALL  Swap(y(i), y(Location))
      END DO
   END SUBROUTINE  Sort



end module mod_sort
