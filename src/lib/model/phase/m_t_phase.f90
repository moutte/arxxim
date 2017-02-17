module M_T_Phase
  use M_Kinds
  implicit none
  !
  private
  !
  public:: T_Phase
  !
  public:: Phase_Zero
  public:: Phase_Index
  public:: Phase_Calc
  !
  type:: T_Phase
  ! describes a phase, including pure phases and mixture or solution phases
    character(len=23):: NamFs
    !! character(len=4) :: Typ !"PURE","AQUEOUS","MIXTURE"
    !! ! "PURE" includes also "DISCRET" 
    !! ! Typ is redundant with the values of iSpc/iMix/iSol
    !! ! -> should not be included in the type declaration
    !! ! pure     == (iSpc/=0, iMix=0,  iSol=0 )
    !! ! MIXTURE  == (iSpc=0,  iMix/=0, iSol=0 )
    !! ! SOLUTION == (iSpc=0,  iMix=0,  iSol/=0)
    integer :: iSpc != index of species in vSpc     -> for pure phase 
    integer :: iMix != index of mixture in vMixFas  -> for mixture phase
    integer :: iSol != index of solution in vSolFas -> for solution phase
    !
    ! variable parameters (intensive, e.g. molar, or specific)
    real(dp):: Grt,VolM3,WeitKg !,H,S,Cp 
    ! molar values at current (T,P), for the phase
    ! Grt=    G/RT
    ! VolM3=  molar volume, M3/Mole
    ! WeitKg= molar weight, Kg/Mole
    !
    real(dp):: MolFs !extensive, mole number of phase in system
    !-> should be "outside" the structure ??
  end type T_Phase
  !
  !type:: T_PhaseDat
  !  real(dp):: Grt,VolM3,WeitKg
  !  real(dp):: Mole
  !end type T_FasData
  !
  type(T_Phase):: Phase_Zero= T_Phase( &
  & "Z",          & !Name,Typ
  & 0,0,0,        & !iSpc,iMix,iSol
  & Zero,One,One, & !Grt,V,WeitKg
  & Zero)           !Mole
  !
contains

!subroutine Phase_Discrete_Init(vDisFas)
!  type(T_Phase),intent(out):: vDisFas(:)
!  !
!end subroutine Phase_Discrete_Init

subroutine Phase_Calc( &
& TdgK,Pbar,vSpc,vMixModel,vMixFas, &
& Fas)
!--
!-- calculate molar properties (weight, volume, Gibbs, ...) of phase Fas
!--
  use M_T_Species, only: T_Species
  use M_T_MixModel,only: T_MixModel
  use M_T_MixPhase,only: T_MixPhase
  use M_T_MixPhase,only: MixPhase_Weit,MixPhase_GibbsRT,MixPhase_Volume
  !
  real(dp),        intent(in)   :: TdgK,Pbar
  type(T_Species), intent(in)   :: vSpc(:)
  type(T_MixModel),intent(in)   :: vMixModel(:)
  type(T_MixPhase),intent(in)   :: vMixFas(:)
  type(T_Phase),   intent(inout):: Fas
  !
  integer :: I
  
  if(Fas%iSpc /= 0) then !----------------------------update pure phases
    !
    I=Fas%iSpc
    Fas%WeitKg= vSpc(I)%WeitKg
    Fas%Grt=    vSpc(I)%G0rt
    Fas%VolM3=  vSpc(I)%V0
    !print *,"NamFs,iSpc,V0=", &
    !& trim(Fas%NamFs)," ",Fas%iSpc," ",vSpc(I)%V0*1.D6
    !
  else if(Fas%iMix /= 0) then !------------------------update mix phases
    !
    I=Fas%iMix !-> index in list of mixture phases
    !
    Fas%WeitKg= MixPhase_Weit( &
    & vSpc,vMixModel(vMixFas(I)%iModel),vMixFas(I))
    !
    Fas%Grt= MixPhase_GibbsRT( &
    & TdgK,Pbar, &
    & vSpc,vMixModel(vMixFas(I)%iModel),vMixFas(I))
    !
    Fas%VolM3= MixPhase_Volume( &
    & vSpc,vMixModel(vMixFas(I)%iModel),vMixFas(I))
    !
  else if (Fas%iSol /= 0) then
    ! todo
  end if
  
  !! select case(Fas%Typ)
  !! 
  !! case("PURE") !,"DISCRET")
  !!   I=Fas%iSpc
  !!   
  !!   Fas%WeitKg= vSpc(I)%WeitKg
  !!   Fas%Grt=    vSpc(I)%G0rt
  !!   Fas%VolM3=  vSpc(I)%V0
  !! 
  !! case("MIXT") 
  !!   I=Fas%iMix !-> index in list of mixture phases
  !!   !
  !!   Fas%WeitKg= MixPhase_Weit( &
  !!   & vSpc,vMixModel(vMixFas(I)%iModel),vMixFas(I))
  !!   !
  !!   Fas%Grt= MixPhase_GibbsRT( &
  !!   & TdgK,Pbar, &
  !!   & vSpc,vMixModel(vMixFas(I)%iModel),vMixFas(I))
  !!   !
  !!   Fas%VolM3= MixPhase_Volume( &
  !!   & vSpc,vMixModel(vMixFas(I)%iModel),vMixFas(I))
  !! 
  !! end select
  
  return
end subroutine Phase_Calc

!! subroutine Phase_Zero(Fas)
!!   !
!!   type(T_Phase),intent(out):: Fas
!!   !
!!   Fas%Name=   "Z"
!!   Fas%Typ=    "PURE"
!!   Fas%iSpc=   0
!!   Fas%iMix=   0
!!   Fas%iSol=   0
!!   Fas%Grt=    Zero
!!   Fas%V=      One
!!   Fas%WeitKg= One
!!   !
!!   Fas%Mole=   Zero
!!   !Fas%H=      Zero
!!   !Fas%S=      Zero
!!   !Fas%Cp=     Zero
!! end subroutine Phase_Zero

integer function Phase_Index(Str,V)
!--
!-- position of phase named Str in V(1:size(V))
!--
  type(T_Phase),intent(in)::V(:)
  character(*), intent(in)::Str
  !
  integer     ::I
  
  Phase_Index=0  !if Str not found -> I=0
  if(size(V)==0) return
  
  I=0
  do
    I=I+1
    if(trim(Str)==trim(V(I)%NamFs)) then
      Phase_Index=I ; exit
    end if
    if(I==size(V)) exit
  end do
  
  return
end function Phase_Index

end module M_T_Phase

