module M_T_KinFas 
!--
!-- data structure for all non-aqueous phases contained in the box:
!-- kinetic phases, equilibrium phases, inert phases, adsoprtion sites
!--
  use M_Kinds
  use M_T_Kinmodel,only: T_KinModel
  implicit none
  !
  private
  !
  public:: T_KinFas
  public:: T_KinFas_Dat
  !! public:: T_KinModel
  !
  public:: KinFas_InitTexture
  public:: KinFas_Zero
  public:: KinFas_Surf_Zero
  !
  public:: KinFas_CalcSurf
  public:: KinFas_CalcSurf_Crunch
  !
  type:: T_KinFas_Dat
  !container for variable parameters of a kinetic phase
    real(dp) ::  &   !
    & SferNumber,&   !nr of grains, fixed when cModel=="D"
    & Radius,    &   !grain radius
    & Surf,      &   !surface / (volume rock+pore)
    & SurfKg         !specific surface, m2/ kg of kin'phase
    real(dp) :: PhiM !volume fraction of mineral IN ROCK (input), then IN BOX (runtime)
    real(dp) :: QsK  !current value of saturation index
    character:: cSat !-> saturation state; cSatur=D(ISSOLU,P(RECIPI,M(INIMAL,I(NERT
  end type T_KinFas_Dat
  !
  type:: T_KinFas
  !container for a kinetic phase (generally minerals) with kinetics and texture data
    character(len=23):: NamKF !name of kinetic phase
    integer:: iFas !index of phase in phase list vFas
    integer:: iKin !index of phase in list of kinetic models, vKinModel
    character:: cMode
    !D-> Density
    !R-> Radius
    !S-> Surface
    !E-> Equilibrium
    !I-> Inert
    !A-> Adsorption site
    real(dp) :: ReacCoef !reactive surface coefficient
    real(dp) :: QsKSeuil !supersaturation threshold (relative ?)
    !! real(dp) :: Surf0    !
    !
    !contains parameters that may vary during evolution
    type(T_KinFas_Dat):: Dat
  end type T_KinFas
  !
contains

subroutine KinFas_Zero(M)
  type(T_KinFas),intent(inout):: M
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
end subroutine KinFas_Zero

subroutine KinFas_Surf_Zero(vFas,M)
!-- initialize surface parameters
!-- given the grain size (M%Dat%Radius),
!-- or given the specific surface (M%Dat%SurfKg),
!-- calculate either M%Dat%SurfKg or M%Dat%Radius
!--
  use M_T_Phase, only: T_Phase 
  !
  type(T_Phase), intent(in)   :: vFas(:)
  type(T_KinFas),intent(inout):: M
  !
  real(dp):: VolKg
  !
  VolKg= vFas(M%iFas)%VolM3 /vFas(M%iFas)%WeitKg
  !
  !Sgrain=4pi.r2 ; Vgrain=4/3pi.r3
  if(M%Dat%Radius>Zero) then
    M%Dat%SurfKg= 3.0D0 *VolKg /M%Dat%Radius
  end if
  !
  if(M%Dat%SurfKg>Zero) then
    M%Dat%Radius= 3.0D0 *VolKg /M%Dat%SurfKg
  end if
  !
end subroutine KinFas_Surf_Zero

subroutine KinFas_InitTexture(  &
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
  use M_Numeric_Const,  only: Pi
  use M_T_Phase,only: T_Phase 
  !
  type(T_Phase), intent(in) :: vFas(:)
  type(T_KinFas),intent(inout):: M
  real(dp),      intent(in) :: Vol        !volume of phase M
  !
  real(dp),      intent(out):: Surf       !surface of volume V of M
  real(dp),      intent(out):: SferNumber !nr.grains in volume V of M
  !
  real(dp):: VGrain,SGrain,Weit
  !
  Weit= Vol /vFas(M%iFas)%VolM3 *vFas(M%iFas)%WeitKg
  !
  if(M%Dat%Radius>Zero) then
    SGrain=       4.0D0*Pi *M%Dat%Radius*M%Dat%Radius
    VGrain=       SGrain   *M%Dat%Radius /3.0D0
    Surf=         SGrain *Vol/VGrain
    M%Dat%SurfKg= Surf /Weit
  end if
  !
  if(M%Dat%SurfKg>Zero) then
    Surf=         M%Dat%SurfKg *Weit
    M%Dat%Radius= Vol /Surf *3.0D0
    SGrain=       4.0D0*Pi *M%Dat%Radius*M%Dat%Radius
    VGrain=       SGrain   *M%Dat%Radius /3.0D0
    !! SferNumber=   Vol/VGrain !Nr.Grains of radius R in volume Vol
  end if
  !
  SferNumber=   Vol/VGrain !Nr.Grains of radius R in volume Vol
  !
end subroutine KinFas_InitTexture

subroutine KinFas_CalcSurf(& !surface of nMol of phase M
!--
!-- obsolete : better use method KinFas_Surf in module M_KinFas_Surf
!--
& vFas,     & !IN:     used for molar volume
& bImplicit,& !IN:
& nMol,     & !IN:     Nr Moles of phase
!! & SurfMinim,& !IN:
& M,        & !IN/OUT: Mineral (update Radius or SferNumber, depending on cModel)
& Surf,     & !OUT:    Surface
& dSRdLnX_, & !OUT:    dSRdLnX_=dSRdLnXm(iMk,iMk)
& dSRdX_)     !OUT:    (for future use ?)
  !
  use M_Numeric_Const,only:PI
  use M_T_Phase,only: T_Phase 
  !
  type(T_Phase), intent(in)   :: vFas(:)
  logical,       intent(in)   :: bImplicit
  real(dp),      intent(in)   :: nMol
  !! real(dp),      intent(in)   :: SurfMinim
  !mod 14/06/2008 12:44 removed SurfMinim
  type(T_KinFas),intent(inout):: M !in=Mole,SferNumber/Radius,VMol, out=Surf
  real(dp),      intent(out)  :: Surf, dSRdLnX_, dSRdX_  !
  !
  real(dp):: Vol,Weit
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
  if (nMol>Zero) then
    Vol=  nMol *vFas(M%iFas)%VolM3
    Weit= nMol *vFas(M%iFas)%WeitKg
    select case(M%cMode)
    case("D")
      ! constant sphere density = growth only
      !(SferNumber=number of grains of mineral M in the box)
      !Surf=GrainNr*SGrain=GrainNr*(4pi*R^2)
      !   = GrainNr*4pi*(3/4pi)^2/3*(Vol/GrainNr)^2/3
      !   = (9*4pi)^1/3 *GrainNr^1/3 *(X.VMol)^2/3 
      Surf=(36*Pi *M%Dat%SferNumber)**(One/3.0D0) *Vol**(2.0D0/3.0D0) !->real surface
      M%Dat%Radius=3.0D0*Vol /Surf !update sphere radius, R=3.Vol/Surf
      !
      if(bImplicit) then
        !Surf=Surf*M%ReacCoef !reactive surface
        dSRdLnX_=(2.0D0/3.0D0)*Surf !dSR/dLnX=X*dSR/dX !->dSRdLnXm(iMn,iMn)
        dSRdX_=  (2.0D0/3.0D0)*Surf /nMol !->dSRdXm(iMn,iMn)
        !else no need to calculate derivative
      end if
    case("R")
      ! constant sphere radius = nucleation only
      Surf= Vol *3.0D0 /M%Dat%Radius !real surface,= MolNumb*VMol*3/R
      M%Dat%SferNumber=Surf/(4*Pi *M%Dat%Radius**2) !update sphere density !!!!!!!!check
      !
      !Surf=Surf*M%ReacCoef !reactive surface
      if(bImplicit) then
        dSRdLnX_= Surf        !->dSRdLnXm(iMn,iMn)
        dSRdX_=   Surf/nMol   !->dSRdXm(iMn,iMn)
        !else no need to calculate derivative
      end if
    end select
    !udpate specific surface (M2/KG)
    M%Dat%SurfKg= Surf /Weit
  end if
  !! if(Surf<SurfMinim) then
  !!   Surf= SurfMinim !SurfMinim= minimal surface per m^3
  !!   dSRdLnX_= Zero !! ?? (2.0D0/3.0D0)*Surf ??
  !!   dSRdX_=   Zero
  !! end if
  !??? needs correction for VBox ???
end subroutine KinFas_CalcSurf

!-----------------------------------------------------------------------
!= subroutine KinFas_CalcSurf_Crunch ===================================
!-----------------------------------------------------------------------
!= OBSOLETE >> use method KinFas_Surf_Crunch in module M_KinFas_Surf
!=
!= derived from method described in CRUNCHFLOW manual (Steefel)
!= currently implemented only in explicit surface mode
!-----------------------------------------------------------------------
subroutine KinFas_CalcSurf_Crunch( &
  vFas,     &   !IN:
  M,        &   !IN: 
  nMol,     &   !IN: Nr Moles of mineral
  VRef,     &   !IN: the volume (rock + pores) containing the nMol, e.g. VBox  
  PhiF,PhiF0,PhiM0,SurfM0, & !IN
  Surf,SurfKg)  !OUT: Surface of nMol, specific surface
  !! dSRdLnX_, & !OUT:    dSRdLnX_=dSRdLnXm(iMk,iMk)
  !! dSRdX_)     !OUT:    (for future use ?)
  use M_T_Phase,only: T_Phase 
  !
  type(T_Phase), intent(in)   :: vFas(:)
  type(T_KinFas),intent(in)   :: M
  real(dp),      intent(in)   :: nMol 
  real(dp),      intent(in)   :: VRef
  real(dp),      intent(in)   :: PhiF,PhiF0,PhiM0,SurfM0
  !
  real(dp),      intent(out)  :: Surf
  real(dp),      intent(out)  :: SurfKg
  !real(dp),     intent(out)  :: dSRdLnX_, dSRdX_  !
  !
  real(dp):: PhiM,Weit
  !
  PhiM= nMol *vFas(M%iFas)%VolM3 /VRef
  Weit= nMol *vFas(M%iFas)%WeitKg
  !
  select case (M%Dat%cSat)
  case("D") !DISSOLU
    Surf= SurfM0 *(PhiF*PhiM/PhiF0/PhiM0)**(2.0D0/3.0D0)
    !in CRUNCH, dissolution law for secondary minerals is
    !Surf= M%Surf0 PhiF*PhiM/PhiF0
  case("P") !PRECIPI
    Surf= SurfM0 *(PhiF/PhiF0)**(2.0D0/3.0D0)
  end select
  !udpate specific surface (M2/KG)
  SurfKg= Surf /Weit
  !
  ! if(Surf<SurfMinim) then
  !   Surf= SurfMinim !SurfMinim= minimal surface per m^3
  !   !! dSRdLnX_= Zero !! ?? (2.0D0/3.0D0)*Surf ??
  !   !! dSRdX_=   Zero
  ! end if
  !??? needs correction for VBox ???
end subroutine KinFas_CalcSurf_Crunch

end module M_T_KinFas

