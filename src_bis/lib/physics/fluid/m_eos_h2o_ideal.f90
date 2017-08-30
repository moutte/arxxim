module M_Eos_H2O_Ideal
!! was M_Ideal_H2O

!-------------------------------------------------------------------------------
! H2O PERFECTGAS : Eos_H2O_Ideal
!-------------------------------------------------------------------------------
! Purpose : Thermodynamic properties for H2O in the ideal gas state
! 
! Reference : SupCrt, equations given by Woolley (1979) 
!             SupCrt92 Subroutine IDEAL [ h2o.f ]
!
!-------------------------------------------------------------------------------
! Arxim Integration : J.Moutte
!-------------------------------------------------------------------------------

use M_Kinds
implicit none
private

real(dp), parameter :: R_jK= 8.314510D0

!// public functionS
public:: Eos_H2O_Ideal ! was CalcH2O_Ideal

contains

subroutine Eos_H2O_Ideal(T,AI,GI,SI,UI,HI,CVI,CPI)
!--
!-- retrieved from SupCrt (routine called ideal),
!-- thermodynamic properties for H2O (AI,GI,SI,UI,HI,CVI,CPI)
!-- in the ideal gas state 
!-- using equations given by Woolley (1979)
!--
  implicit none
  real(dp),intent(in) ::T
  real(dp),intent(out)::AI,GI,SI,UI,HI,CVI,CPI
  real(dp)::&
    C(18), &
    TT,TL,EMULT
  integer::I
  !
  C=(/ &
  &  0.197302710180D02,  0.209662681977D2,  -0.483429455355D0,  &
  &  0.605743189245D01,  0.225602388500D2,  -0.987532442000D1,  &
  & -0.431355385130D01,  0.458155781000D0,  -0.477549018830D-1, &
  &  0.412384606330D-2, -0.279290528520D-3,  0.144816952610D-4, &
  & -0.564736587480D-6,  0.162004460000D-7, -0.330382279600D-9, &
  &  0.451916067368D-11,-0.370734122710D-13, 0.137546068238D-15 /)
  !
  TT= T/100._dp
  TL= log(TT)
  GI= -(C(1)/TT + C(2)) *TL
  HI=  (C(2) + C(1)*(One - TL)/TT)
  CPI=  C(2) - C(1)/TT
  do I=3,18
    EMULT= TT**(I-6)
    GI=  GI  - C(I) *EMULT
    HI=  HI  + C(I) *(I-6) *EMULT
    CPI= CPI + C(I) *(I-6) *(I-5) *EMULT
  enddo 
  AI= GI - One
  UI= HI - One
  CVI= CPI - One
  SI= UI - AI
  !
  return
end subroutine Eos_H2O_Ideal

end module M_Eos_H2O_Ideal
