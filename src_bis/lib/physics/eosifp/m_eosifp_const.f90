module M_T_EosIfp_Const

implicit none
private

!// General Cubic Laws Parameters 
type:: T_EosIfp_Const

  !// Common Parameters PrModel and SrkModel
  real(kind=8) :: XAF(0:3,4)    ! Common Cubic Parameters PR and SRK

  !// Common Parameters PrModel
  real(kind=8) :: OMGAPR        ! Omega A coeff PR
  real(kind=8) :: OMGBPR        ! Omega B coeff PR
  real(kind=8) :: DLPR1         ! DL PR 1
  real(kind=8) :: DLPR2         ! DL PR 2

  !// Common Parameters SrkModel
  real(kind=8) :: OMGARK        ! Omega A coeff SRK
  real(kind=8) :: OMGBRK        ! Omega B coeff SRK
  real(kind=8) :: DLRK1         ! DL SRK 1
  real(kind=8) :: DLRK2         ! DL SRK 2

end type T_EosIfp_Const

!// ============================
public :: T_EosIfp_Const 
public :: EosIfp_Const_init

contains

subroutine EosIfp_Const_init( T )

  implicit none
  type(T_EosIfp_Const), intent(inout) :: T

  real(kind=8), parameter :: RAC2 = 1.4142136562d0

  !------- data THERMO

  T% XAF (0,:) = 1.d0 * (/ .37464 , 1.54226, -.26992 , .0 /)
  T% XAF (1,:) = 1.d0 * (/ .379642, 1.48503, -.164423, .016666 /)
  T% XAF (2,:) = 1.d0 * (/ .480   , 1.574  , -.176   , .0 /)
  T% XAF (3,:) = 1.d0 * (/ .48508 , 1.55171, -.15613 , .0 /)

  !------ data PENG ROBINSON

  T% OMGAPR = 0.457235d0
  T% OMGBPR = 0.077796d0      
  T% DLPR1 =  1. - RAC2
  T% DLPR2 =  1. + RAC2    

  ! ===== doNNEES COMMUNES SRK

  T% OMGARK = 0.42747d0
  T% OMGBRK = 0.08664d0
  T% DLRK1 =  0.
  T% DLRK2 =  1. 

end subroutine EosIfp_Const_init

end module M_T_EosIfp_Const
