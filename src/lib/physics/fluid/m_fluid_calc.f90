MODULE M_Fluid_Calc

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

  USE M_Kinds

  !// Import Modules Fluid H2O
  USE M_Eos_H2O_Rho_Psat
  USE M_Eos_H2O_Ideal
  USE M_Eos_H2O_Holland_Powell
  USE M_Eos_H2O_Haar
  USE M_Eos_H2O_Haar_Ghiorso
  
  !// Import Modules Fluid CO2
  USE M_Eos_CO2_Holland_Powell
  USE M_Eos_KerJac
  USE M_Mixmodel_Duan
  
  IMPLICIT NONE

  PRIVATE

  !// EXPORT FUNCTIONS FROM : M_Eos_H2O_Rho_Psat
  PUBLIC:: Eos_H2O_rho
  PUBLIC:: Eos_H2O_psat

  !// EXPORT FUNCTIONS FROM : M_Eos_H2O_Ideal
  PUBLIC:: Eos_H2O_Ideal

  !// EXPORT FUNCTIONS FROM : M_Eos_H2O_Holland_Powell
  PUBLIC:: Eos_H2O_HolPow91
  PUBLIC:: Eos_H2O_HolPow

  !// EXPORT FUNCTIONS FROM : M_Supcrt_H2O
  ! PUBLIC:: CalcGH2O_Supcrt

  !// EXPORT FUNCTIONS FROM : M_Eos_H2O_Haar
  PUBLIC:: Eos_H2O_Haar

  !// EXPORT FUNCTIONS FROM : M_Eos_H2O_Haar_Ghiorso
  PUBLIC:: Eos_H2O_Haar_Ghiorso

  !// EXPORT FUNCTIONS FROM : M_Eos_CO2_HolPow91
  PUBLIC:: Eos_CO2_HolPow91

  !// PUBLIC FUNCTION WRAPPER
  PUBLIC:: CalcFluid

  CONTAINS

SUBROUTINE CalcFluid( &
& sFluid,sModel, &
& TdgK,Pbar,     &
& Grt,H,S,V_m3,  &
& Ok)
  CHARACTER(LEN=*),INTENT(IN) :: sFluid,sModel
  REAL(dp),        INTENT(IN) :: TdgK,Pbar
  !
  REAL(dp),        INTENT(OUT):: Grt,H,S,V_m3
  LOGICAL,         INTENT(OUT):: Ok

  Ok= .TRUE.

  Grt= Zero
  H=   Zero
  S=   Zero
  V_m3=Zero
  !
  SELECT CASE(TRIM(sFluid))
  !
  CASE("H2O")

    SELECT CASE(TRIM(sModel))
    CASE("HAAR.GHIORSO")
      CALL Eos_H2O_Haar_Ghiorso(TdgK,Pbar,Grt,H,S,V_m3)
    !CASE("SUPCRT")
    !  CALL CalcGH2O_Supcrt(TdgK,Pbar,Grt,H,S,V_m3)
    CASE("HAAR")
      CALL Eos_H2O_Haar  (TdgK,Pbar,Grt,V_m3)
    CASE("HP91")
      CALL Eos_H2O_HolPow91  (TdgK,Pbar,Grt,V_m3)
    CASE("HOLPOW")
      CALL Eos_H2O_HolPow(TdgK,Pbar,Grt)
    !CASE("HAAR2")
      !CALL WHAAR2(Pbar,TdgK,Grt,V_m3)
    CASE DEFAULT
      !! err_msg= TRIM(sFluid) // "=unknown model for H2O !!"
      Ok= .FALSE.
    END SELECT

  CASE("CO2")

    SELECT CASE(TRIM(sModel))
    CASE("HP91")
      CALL Eos_CO2_HolPow91(TdgK,Pbar,Grt,V_m3)
    END SELECT

  END SELECT

ENDSUBROUTINE CalcFluid

ENDMODULE M_Fluid_Calc
