module M_KinFas_Surf
!--
!-- compute mineral surface and derivatives for kinetic rate laws
!--
  use M_Kinds
  implicit none
  !
  private
  !
  public:: KinFas_Surf_Init
  public:: KinFas_Surf_Update
  public:: KinFas_Surf_Calc
  public:: KinFas_Surf_Crunch
  !
contains

subroutine KinFas_Surf_Init( &
!-- initialize surface -------------------------------------------------
!-- given the specific surface (M%Dat%SurfKg),
!--.compute surface for a given volume
!-----------------------------------------------------------------------
& vFas,      & !
& M,         & !
& Vol,       & !in:  volume of phase M
& Surf)        !out: surface of volume V of M
  use M_Numeric_Const,  only: Pi
  use M_T_KinFas,only: T_KinFas
  use M_T_Phase, only: T_Phase 
  !
  type(T_Phase), intent(in):: vFas(:)
  type(T_KinFas),intent(in):: M
  real(dp),      intent(in):: Vol        !volume of phase M
  !
  real(dp),      intent(out):: Surf       !surface of volume V of M
  !
  Surf= M%Dat%SurfKg *Vol *vFas(M%iFas)%WeitKg /vFas(M%iFas)%VolM3
  !
end subroutine KinFas_Surf_Init

subroutine KinFas_Surf_Update( &
!-- update surface parameters-------------------------------------------
!-- given the mole number and surface of phase,
!-- compute specific surface and equivalent radius
!-----------------------------------------------------------------------
& vFas,      & !
& nMol,      & !in: nr mole of phase M
& Surf,      & !in: surface of phase M
& M)           !inout
  use M_Numeric_Const,  only: Pi
  use M_T_KinFas,only: T_KinFas
  use M_T_Phase, only: T_Phase 
  !
  type(T_Phase), intent(in) :: vFas(:)
  real(dp),      intent(in) :: nMol, Surf
  type(T_KinFas),intent(inout):: M
  !
  M%Dat%SurfKg= Surf /nMol/vFas(M%iFas)%WeitKg
  !
  M%Dat%Radius= 3.0D0 *nMol*vFas(M%iFas)%VolM3 /Surf 
  !
end subroutine KinFas_Surf_Update

subroutine KinFas_Surf_Calc(& !
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
& dSRdX_)     !OUT: (for future use ?)
  !
  logical,  intent(in) :: bImplicit
  character,intent(in) :: cModel
  real(dp), intent(in) :: nMol0, nMol, Surf0
  real(dp), intent(out):: Surf, dSRdLnX_, dSRdX_  !
  !
  Surf=     Zero
  dSRdLnX_= Zero
  dSRdX_=   Zero
  !
  !write(73,'(2G15.6)') Surf0,nMol0
  select case(cModel)
  case("D")
    !--------------------- D (constant sphere Density) - pure growth ---
    Surf= Surf0 *(nMol /nMol0)**(2.0D0/3.0D0) !-> "geometric" surface
    !Surf= Surf*M%ReacCoef !reactive surface
    if(bImplicit) then
      dSRdLnX_= 2.0D0 /3.0D0 *Surf !dSR/dLnX=X*dSR/dX !->dSRdLnXm(iMn,iMn)
      dSRdX_=   2.0D0 /3.0D0 *Surf /nMol !->dSRdXm(iMn,iMn)
      !else no need to calculate derivative
    end if
  case("R")
    !------------------ R (constant sphere Radius) - pure nucleation ---
    Surf=  Surf0 *nMol /nMol0
    !Surf= Surf*M%ReacCoef !reactive surface
    if(bImplicit) then
      dSRdLnX_= Surf           !->dSRdLnXm(iMn,iMn)
      dSRdX_=   Surf0 /nMol0   !->dSRdXm(iMn,iMn)
      !else no need to calculate derivative
    end if
  end select
  !
end subroutine KinFas_Surf_Calc

subroutine KinFas_Surf_Crunch( &
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
!! dSRdX_)      !OUT:    (for future use ?)
  !
  character,intent(in) :: cSat
  real(dp), intent(in) :: nMol,nMol0,Surf0
  real(dp), intent(in) :: PhiF,PhiF0
  !
  real(dp), intent(out):: Surf
  !real(dp),     intent(out)  :: dSRdLnX_, dSRdX_  !
  !
  select case(cSat)
  case("D") !"DISSOLU"
    Surf= Surf0 *(nMol/nMol0)**(2.0D0/3.0D0) *(PhiF/PhiF0)**(2.0D0/3.0D0)
    ! or Surf= (nMol*PhiF)**(2.0D0/3.0D0) *Surf0 / (nMol0*PhiF0)**(2.0D0/3.0D0)
    !in CRUNCH, dissolution law for secondary minerals is
    !Surf= M%Surf0 PhiF*PhiM/PhiF0
  case("P") !"PRECIPI"
    Surf= Surf0 *(PhiF/PhiF0)**(2.0D0/3.0D0)
  end select
  !
end subroutine KinFas_Surf_Crunch

end module M_KinFas_Surf

