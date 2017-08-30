module M_Fluid_Calc

!-----------------------------------------------------------------------
! Purpose : General Interface for Computation of Pure Fluid Properties
!
! > Pure H2O (Water , Steam)
! > Pure CO2 (Liquid, Vapor)
!
!-----------------------------------------------------------------------
!
! calculations of HSV Properties of Fluid : GH2O_Haar, GH2O_HolPw, ...
! partly adapted from SupCrt92
! and from CdeCapitani F77 sources
! (http://titan.minpet.unibas.ch/minpet/theriak/)
!-----------------------------------------------------------------------

  use M_Kinds

  !// Import Modules Fluid H2O
  use M_Eos_H2O_Rho_Psat
  use M_Eos_H2O_Ideal
  use M_Eos_H2O_Holland_Powell
  use M_Eos_H2O_Haar
  use M_Eos_H2O_Haar_Ghiorso
  
  !// Import Modules Fluid CO2
  use M_Eos_CO2_Holland_Powell
  use M_Eos_KerJac
  use M_Mixmodel_Duan
  
  implicit none

  private

  !// EXPORT functionS FROM : M_Eos_H2O_Rho_Psat
  public:: Eos_H2O_rho
  public:: Eos_H2O_psat

  !// EXPORT functionS FROM : M_Eos_H2O_Ideal
  public:: Eos_H2O_Ideal

  !// EXPORT functionS FROM : M_Eos_H2O_Holland_Powell
  public:: Eos_H2O_HolPow91
  public:: Eos_H2O_HolPow

  !// EXPORT functionS FROM : M_Supcrt_H2O
  ! public:: CalcGH2O_Supcrt

  !// EXPORT functionS FROM : M_Eos_H2O_Haar
  public:: Eos_H2O_Haar

  !// EXPORT functionS FROM : M_Eos_H2O_Haar_Ghiorso
  public:: Eos_H2O_Haar_Ghiorso

  !// EXPORT functionS FROM : M_Eos_CO2_HolPow91
  public:: Eos_CO2_HolPow91

  !// public function WRAPPER
  public:: CalcFluid

  contains

subroutine CalcFluid( &
& sFluid,sModel, &
& TdgK,Pbar,     &
& Grt,H,S,V_m3,  &
& Ok)
  character(len=*),intent(in) :: sFluid,sModel
  real(dp),        intent(in) :: TdgK,Pbar
  !
  real(dp),        intent(out):: Grt,H,S,V_m3
  logical,         intent(out):: Ok

  Ok= .true.

  Grt= Zero
  H=   Zero
  S=   Zero
  V_m3=Zero
  !
  select case(trim(sFluid))
  !
  case("H2O")

    select case(trim(sModel))
    case("HAAR.GHIORSO")
      call Eos_H2O_Haar_Ghiorso(TdgK,Pbar,Grt,H,S,V_m3)
    !case("SUPCRT")
    !  call CalcGH2O_Supcrt(TdgK,Pbar,Grt,H,S,V_m3)
    case("HAAR")
      call Eos_H2O_Haar  (TdgK,Pbar,Grt,V_m3)
    case("HP91")
      call Eos_H2O_HolPow91  (TdgK,Pbar,Grt,V_m3)
    case("HOLPOW")
      call Eos_H2O_HolPow(TdgK,Pbar,Grt)
    !case("HAAR2")
      !call WHAAR2(Pbar,TdgK,Grt,V_m3)
    case default
      !! err_msg= trim(sFluid) // "=unknown model for H2O !!"
      Ok= .false.
    end select

  case("CO2")

    select case(trim(sModel))
    case("HP91")
      call Eos_CO2_HolPow91(TdgK,Pbar,Grt,V_m3)
    end select

  end select

end subroutine CalcFluid

end module M_Fluid_Calc
