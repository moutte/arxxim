MODULE M_KinFas_Surf
!--
!-- compute mineral surface and derivatives for kinetic rate laws
!--
  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: KinFas_Surf_Init
  PUBLIC:: KinFas_Surf_Update
  PUBLIC:: KinFas_Surf_Calc
  PUBLIC:: KinFas_Surf_Crunch
  !
CONTAINS

SUBROUTINE KinFas_Surf_Init( &
!-- initialize surface -------------------------------------------------
!-- given the specific surface (M%Dat%SurfKg),
!--.compute surface for a given volume
!-----------------------------------------------------------------------
& vFas,      & !
& M,         & !
& Vol,       & !in:  volume of phase M
& Surf)        !out: surface of volume V of M
  USE M_Numeric_Const,  ONLY: Pi
  USE M_T_KinFas,ONLY: T_KinFas
  USE M_T_Phase, ONLY: T_Phase 
  !
  TYPE(T_Phase), INTENT(IN):: vFas(:)
  TYPE(T_KinFas),INTENT(IN):: M
  REAL(dp),      INTENT(IN):: Vol        !volume of phase M
  !
  REAL(dp),      INTENT(OUT):: Surf       !surface of volume V of M
  !
  Surf= M%Dat%SurfKg *Vol *vFas(M%iFas)%WeitKg /vFas(M%iFas)%VolM3
  !
ENDSUBROUTINE KinFas_Surf_Init

SUBROUTINE KinFas_Surf_Update( &
!-- update surface parameters-------------------------------------------
!-- given the mole number and surface of phase,
!-- compute specific surface and equivalent radius
!-----------------------------------------------------------------------
& vFas,      & !
& nMol,      & !in: nr mole of phase M
& Surf,      & !in: surface of phase M
& M)           !inout
  USE M_Numeric_Const,  ONLY: Pi
  USE M_T_KinFas,ONLY: T_KinFas
  USE M_T_Phase, ONLY: T_Phase 
  !
  TYPE(T_Phase), INTENT(IN) :: vFas(:)
  REAL(dp),      INTENT(IN) :: nMol, Surf
  TYPE(T_KinFas),INTENT(INOUT):: M
  !
  M%Dat%SurfKg= Surf /nMol/vFas(M%iFas)%WeitKg
  !
  M%Dat%Radius= 3.0D0 *nMol*vFas(M%iFas)%VolM3 /Surf 
  !
ENDSUBROUTINE KinFas_Surf_Update

SUBROUTINE KinFas_Surf_Calc(& !
!--
!-- compute surface of nMol of phase M
!--
& bImplicit,& !IN
& cModel,   & !IN
& nMol,     & !IN : Nr Moles of phase
& nMol0,    & !IN : Nr Moles of phase at init' state
& Surf0,    & !IN : Surface of phase at init' state
& Surf,     & !OUT: Surface
& dSRdLnX_, & !OUT: dSRdLnX_=dSRdLnXm(iMk,iMk)
& dSRdX_)     !OUT: (for future USE ?)
  !
  LOGICAL,  INTENT(IN) :: bImplicit
  CHARACTER,INTENT(IN) :: cModel
  REAL(dp), INTENT(IN) :: nMol0, nMol, Surf0
  REAL(dp), INTENT(OUT):: Surf, dSRdLnX_, dSRdX_  !
  !
  Surf=     Zero
  dSRdLnX_= Zero
  dSRdX_=   Zero
  !
  !WRITE(73,'(2G15.6)') Surf0,nMol0
  SELECT CASE(cModel)
  CASE("D")
    !--------------------- D (constant sphere Density) - pure growth ---
    Surf= Surf0 *(nMol /nMol0)**(2.0D0/3.0D0) !-> "geometric" surface
    !Surf= Surf*M%ReacCoef !reactive surface
    IF(bImplicit) THEN
      dSRdLnX_= 2.0D0 /3.0D0 *Surf !dSR/dLnX=X*dSR/dX !->dSRdLnXm(iMn,iMn)
      dSRdX_=   2.0D0 /3.0D0 *Surf /nMol !->dSRdXm(iMn,iMn)
      !ELSE no need to calculate derivative
    ENDIF
  CASE("R")
    !------------------ R (constant sphere Radius) - pure nucleation ---
    Surf=  Surf0 *nMol /nMol0
    !Surf= Surf*M%ReacCoef !reactive surface
    IF(bImplicit) THEN
      dSRdLnX_= Surf           !->dSRdLnXm(iMn,iMn)
      dSRdX_=   Surf0 /nMol0   !->dSRdXm(iMn,iMn)
      !ELSE no need to calculate derivative
    ENDIF
  END SELECT
  !
ENDSUBROUTINE KinFas_Surf_Calc

SUBROUTINE KinFas_Surf_Crunch( &
!--
!-- same surface / volume relation as implemented in CRUNCH (Steefel)
!-- Surf_m / Surf_m0- (Vol_m /Vol_m0)^(2/3) * (PhiF/PhiF0)^(2/3)
!-- makes Surf_m tend to Zero when porosity tends to Zero
!-- currently implemented only in explicit surface mode
!--
& cSat,       &   !IN: 
& nMol,       &   !IN: Nr Moles of phase
& nMol0,      &   !IN: Nr Moles of phase at given time
& Surf0,      &   !IN: Surface at given time
& PhiF,PhiF0, &   !IN: porosity / porosity at given time
& Surf)           !OUT: Surface of nMol
!! dSRdLnX_,&   !OUT:    dSRdLnX_=dSRdLnXm(iMk,iMk)
!! dSRdX_)      !OUT:    (for future USE ?)
  !
  CHARACTER,INTENT(IN) :: cSat
  REAL(dp), INTENT(IN) :: nMol,nMol0,Surf0
  REAL(dp), INTENT(IN) :: PhiF,PhiF0
  !
  REAL(dp), INTENT(OUT):: Surf
  !REAL(dp),     INTENT(OUT)  :: dSRdLnX_, dSRdX_  !
  !
  SELECT CASE(cSat)
  CASE("D") !"DISSOLU"
    Surf= Surf0 *(nMol/nMol0)**(2.0D0/3.0D0) *(PhiF/PhiF0)**(2.0D0/3.0D0)
    ! or Surf= (nMol*PhiF)**(2.0D0/3.0D0) *Surf0 / (nMol0*PhiF0)**(2.0D0/3.0D0)
    !in CRUNCH, dissolution law for secondary minerals is
    !Surf= M%Surf0 PhiF*PhiM/PhiF0
  CASE("P") !"PRECIPI"
    Surf= Surf0 *(PhiF/PhiF0)**(2.0D0/3.0D0)
  END SELECT
  !
ENDSUBROUTINE KinFas_Surf_Crunch

ENDMODULE M_KinFas_Surf

