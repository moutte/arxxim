MODULE M_T_KinFas 
!--
!-- data structure for all non-aqueous phases contained in the box:
!-- kinetic phases, equilibrium phases, inert phases, adsoprtion sites
!--
  USE M_Kinds
  USE M_T_Kinmodel,ONLY: T_KinModel
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: T_KinFas
  PUBLIC:: T_KinFas_Dat
  !! PUBLIC:: T_KinModel
  !
  PUBLIC:: KinFas_InitTexture
  PUBLIC:: KinFas_Zero
  PUBLIC:: KinFas_Surf_Zero
  !
  PUBLIC:: KinFas_CalcSurf
  PUBLIC:: KinFas_CalcSurf_Crunch
  !
  TYPE:: T_KinFas_Dat
  !container for variable parameters of a kinetic phase
    REAL(dp) ::  & 
    & SferNumber,&   !nr of grains, fixed when cModel=="D"
    & Radius,    &   !grain radius
    & Surf,      &   !surface / (volume rock+pore)
    & SurfKg         !specific surface, m2/ kg of kin'phase
    REAL(dp) :: PhiM !volume fraction of mineral IN ROCK (input), THEN IN BOX (runtime)
    REAL(dp) :: QsK  !current value of saturation index
    CHARACTER:: cSat !-> saturation state; cSatur=D(ISSOLU,P(RECIPI,M(INIMAL,I(NERT
  ENDTYPE T_KinFas_Dat
  !
  TYPE:: T_KinFas
  !container for a kinetic phase (generally minerals) with kinetics and texture data
    CHARACTER(LEN=23):: NamKF !name of kinetic phase
    INTEGER:: iFas !index of phase in phase list vFas
    INTEGER:: iKin !index of phase in list of kinetic models, vKinModel
    CHARACTER:: cMode
    !D-> Density
    !R-> Radius
    !S-> Surface
    !E-> Equilibrium
    !I-> Inert
    !A-> Adsorption site
    REAL(dp) :: ReacCoef !reactive surface coefficient
    REAL(dp) :: QsKSeuil !supersaturation threshold (relative ?)
    !! REAL(dp) :: Surf0    !
    !
    !CONTAINS PARAMETERs that may vary during evolution
    TYPE(T_KinFas_Dat):: Dat
  ENDTYPE T_KinFas
  !
CONTAINS

SUBROUTINE KinFas_Zero(M)
  TYPE(T_KinFas),INTENT(INOUT):: M
  !
  M%NamKF=    "X"
  M%iFas=     0
  M%iKin=     0
  M%ReacCoef= One
  M%QsKSeuil= One
  M%cMode=    "D"
  !
  M%Dat%cSat=     "I" !INERT
  M%Dat%SferNumber= 1.0D15
  M%Dat%Radius=   1.0D-6
  M%Dat%SurfKg=   Zero
  M%Dat%PhiM=     Zero
  M%Dat%QsK=      One
  !
ENDSUBROUTINE KinFas_Zero

SUBROUTINE KinFas_Surf_Zero(vFas,M)
!-- initialize surface parameters
!-- given the grain size (M%Dat%Radius),
!-- or given the specific surface (M%Dat%SurfKg),
!-- calculate either M%Dat%SurfKg or M%Dat%Radius
!--
  USE M_T_Phase, ONLY: T_Phase 
  !
  TYPE(T_Phase), INTENT(IN)   :: vFas(:)
  TYPE(T_KinFas),INTENT(INOUT):: M
  !
  REAL(dp):: VolKg
  !
  VolKg= vFas(M%iFas)%VolM3 /vFas(M%iFas)%WeitKg
  !
  !Sgrain=4pi.r2 ; Vgrain=4/3pi.r3
  IF(M%Dat%Radius>Zero) THEN
    M%Dat%SurfKg= 3.0D0 *VolKg /M%Dat%Radius
  ENDIF
  !
  IF(M%Dat%SurfKg>Zero) THEN
    M%Dat%Radius= 3.0D0 *VolKg /M%Dat%SurfKg
  ENDIF
  !
ENDSUBROUTINE KinFas_Surf_Zero

SUBROUTINE KinFas_InitTexture(  &
!-- initialize surface parameters
!-- given the grain size (M%Dat%Radius),
!-- or given the specific surface (M%Dat%SurfKg),
!-- calculate SferNumber and Surface for a volume Vol of phase
!--
& vFas,      & !
& M,         & !
& Vol,       & !volume of phase M
& Surf,      & !surface of volume V of M
& SferNumber)  !nr.grains for a volume V of M 
  USE M_Numeric_Const,  ONLY: Pi
  USE M_T_Phase,ONLY: T_Phase 
  !
  TYPE(T_Phase), INTENT(IN) :: vFas(:)
  TYPE(T_KinFas),INTENT(INOUT):: M
  REAL(dp),      INTENT(IN) :: Vol        !volume of phase M
  !
  REAL(dp),      INTENT(OUT):: Surf       !surface of volume V of M
  REAL(dp),      INTENT(OUT):: SferNumber !nr.grains in volume V of M
  !
  REAL(dp):: VGrain,SGrain,Weit
  !
  Weit= Vol /vFas(M%iFas)%VolM3 *vFas(M%iFas)%WeitKg
  !
  IF(M%Dat%Radius>Zero) THEN
    SGrain=       4.0D0*Pi *M%Dat%Radius*M%Dat%Radius
    VGrain=       SGrain   *M%Dat%Radius /3.0D0
    Surf=         SGrain *Vol/VGrain
    M%Dat%SurfKg= Surf /Weit
  ENDIF
  !
  IF(M%Dat%SurfKg>Zero) THEN
    Surf=         M%Dat%SurfKg *Weit
    M%Dat%Radius= Vol /Surf *3.0D0
    SGrain=       4.0D0*Pi *M%Dat%Radius*M%Dat%Radius
    VGrain=       SGrain   *M%Dat%Radius /3.0D0
    !! SferNumber=   Vol/VGrain !Nr.Grains of radius R in volume Vol
  ENDIF
  !
  SferNumber=   Vol/VGrain !Nr.Grains of radius R in volume Vol
  !
ENDSUBROUTINE KinFas_InitTexture

SUBROUTINE KinFas_CalcSurf(& !surface of nMol of phase M
!--
!-- obsolete : better use method KinFas_Surf in module M_KinFas_Surf
!--
& vFas,     & !IN:     used for molar volume
& bImplicit,& !IN:
& nMol,     & !IN:     Nr Moles of phase
!! & SurfMinim,& !IN:
& M,        & !IN/OUT: Mineral (update Radius or SferNumber, depENDing on cModel)
& Surf,     & !OUT:    Surface
& dSRdLnX_, & !OUT:    dSRdLnX_=dSRdLnXm(iMk,iMk)
& dSRdX_)     !OUT:    (for future USE ?)
  !
  USE M_Numeric_Const,ONLY:PI
  USE M_T_Phase,ONLY: T_Phase 
  !
  TYPE(T_Phase), INTENT(IN)   :: vFas(:)
  LOGICAL,       INTENT(IN)   :: bImplicit
  REAL(dp),      INTENT(IN)   :: nMol
  !! REAL(dp),      INTENT(IN)   :: SurfMinim
  !mod 14/06/2008 12:44 removed SurfMinim
  TYPE(T_KinFas),INTENT(INOUT):: M !in=Mole,SferNumber/Radius,VMol, out=Surf
  REAL(dp),      INTENT(OUT)  :: Surf, dSRdLnX_, dSRdX_  !
  !
  REAL(dp):: Vol,Weit
  !
  !Vol=    R3 *4pi/3
  !Surf=   R2 *4pi
  !Surf^3= Vol^2 * 36pi
  !
  !Vol=nMol*VMol
  !Vol= SferNumber*SferVol
  !Surf=SferNumber*SferSurf
  !Surf/Vol= SferSurf/SferVol= 3/SferRad
  Surf=     Zero
  dSRdLnX_= Zero !! ?? (2.0D0/3.0D0)*Surf ??
  dSRdX_=   Zero
  IF (nMol>Zero) THEN
    Vol=  nMol *vFas(M%iFas)%VolM3
    Weit= nMol *vFas(M%iFas)%WeitKg
    SELECT CASE(M%cMode)
    CASE("D")
      ! constant sphere density = growth only
      !(SferNumber=number of grains of mineral M in the box)
      !Surf=GrainNr*SGrain=GrainNr*(4pi*R^2)
      !   = GrainNr*4pi*(3/4pi)^2/3*(Vol/GrainNr)^2/3
      !   = (9*4pi)^1/3 *GrainNr^1/3 *(X.VMol)^2/3 
      Surf=(36*Pi *M%Dat%SferNumber)**(One/3.0D0) *Vol**(2.0D0/3.0D0) !->REAL surface
      M%Dat%Radius=3.0D0*Vol /Surf !update sphere radius, R=3.Vol/Surf
      !
      IF(bImplicit) THEN
        !Surf=Surf*M%ReacCoef !reactive surface
        dSRdLnX_=(2.0D0/3.0D0)*Surf !dSR/dLnX=X*dSR/dX !->dSRdLnXm(iMn,iMn)
        dSRdX_=  (2.0D0/3.0D0)*Surf /nMol !->dSRdXm(iMn,iMn)
        !ELSE no need to calculate derivative
      ENDIF
    CASE("R")
      ! constant sphere radius = nucleation only
      Surf= Vol *3.0D0 /M%Dat%Radius !REAL surface,= MolNumb*VMol*3/R
      M%Dat%SferNumber=Surf/(4*Pi *M%Dat%Radius**2) !update sphere density !!!!!!!!check
      !
      !Surf=Surf*M%ReacCoef !reactive surface
      IF(bImplicit) THEN
        dSRdLnX_= Surf        !->dSRdLnXm(iMn,iMn)
        dSRdX_=   Surf/nMol   !->dSRdXm(iMn,iMn)
        !ELSE no need to calculate derivative
      ENDIF
    END SELECT
    !udpate specific surface (M2/KG)
    M%Dat%SurfKg= Surf /Weit
  ENDIF
  !! IF(Surf<SurfMinim) THEN
  !!   Surf= SurfMinim !SurfMinim= minimal surface per m^3
  !!   dSRdLnX_= Zero !! ?? (2.0D0/3.0D0)*Surf ??
  !!   dSRdX_=   Zero
  !! ENDIF
  !??? needs correction for VBox ???
ENDSUBROUTINE KinFas_CalcSurf

!-----------------------------------------------------------------------
!= SUBROUTINE KinFas_CalcSurf_Crunch ===================================
!-----------------------------------------------------------------------
!= OBSOLETE >> use method KinFas_Surf_Crunch in module M_KinFas_Surf
!=
!= derived from method described in CRUNCHFLOW manual (Steefel)
!= currently implemented only in explicit surface mode
!-----------------------------------------------------------------------
SUBROUTINE KinFas_CalcSurf_Crunch( &
  vFas,     &   !IN:
  M,        &   !IN: 
  nMol,     &   !IN: Nr Moles of mineral
  VRef,     &   !IN: the volume (rock + pores) containing the nMol, e.g. VBox  
  PhiF,PhiF0,PhiM0,SurfM0, & !IN
  Surf,SurfKg)  !OUT: Surface of nMol, specIFic surface
  !! dSRdLnX_, & !OUT:    dSRdLnX_=dSRdLnXm(iMk,iMk)
  !! dSRdX_)     !OUT:    (for future USE ?)
  USE M_T_Phase,ONLY: T_Phase 
  !
  TYPE(T_Phase), INTENT(IN)   :: vFas(:)
  TYPE(T_KinFas),INTENT(IN)   :: M
  REAL(dp),      INTENT(IN)   :: nMol 
  REAL(dp),      INTENT(IN)   :: VRef
  REAL(dp),      INTENT(IN)   :: PhiF,PhiF0,PhiM0,SurfM0
  !
  REAL(dp),      INTENT(OUT)  :: Surf
  REAL(dp),      INTENT(OUT)  :: SurfKg
  !REAL(dp),     INTENT(OUT)  :: dSRdLnX_, dSRdX_  !
  !
  REAL(dp):: PhiM,Weit
  !
  PhiM= nMol *vFas(M%iFas)%VolM3 /VRef
  Weit= nMol *vFas(M%iFas)%WeitKg
  !
  SELECT CASE (M%Dat%cSat)
  CASE("D") !DISSOLU
    Surf= SurfM0 *(PhiF*PhiM/PhiF0/PhiM0)**(2.0D0/3.0D0)
    !in CRUNCH, dissolution law for secondary minerals is
    !Surf= M%Surf0 PhiF*PhiM/PhiF0
  CASE("P") !PRECIPI
    Surf= SurfM0 *(PhiF/PhiF0)**(2.0D0/3.0D0)
  END SELECT
  !udpate specific surface (M2/KG)
  SurfKg= Surf /Weit
  !
  ! IF(Surf<SurfMinim) THEN
  !   Surf= SurfMinim !SurfMinim= minimal surface per m^3
  !   !! dSRdLnX_= Zero !! ?? (2.0D0/3.0D0)*Surf ??
  !   !! dSRdX_=   Zero
  ! ENDIF
  !??? needs correction for VBox ???
ENDSUBROUTINE KinFas_CalcSurf_Crunch

ENDMODULE M_T_KinFas

