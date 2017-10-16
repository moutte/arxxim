!.parameters for size of real types
module M_Kinds !kind defining the size of variables
  integer, parameter :: dp= KIND(1.0D0) 
  integer, parameter :: sp= KIND(1.0)
  real(dp),parameter :: Zero=0.0_dp, One=1.0_dp, Two=2.0_dp
end module M_Kinds

! To declare single precision, use
!   INTEGER, PARAMETER:: rk = SELECTED_REAL_KIND(6,37)
! To declare a double precision, use
!   INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15,307)
! To declare a Quadruple precision (not very often used):
!   INTEGER, PARAMETER :: rk = selected_real_kind(33, 4931)
!
!   integer, parameter   :: i4=SELECTED_INT_KIND(4)
!   integer, parameter   :: i8=SELECTED_INT_KIND(8)
!
