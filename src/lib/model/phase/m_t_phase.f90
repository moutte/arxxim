MODULE M_T_Phase
  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: T_Phase
  !
  PUBLIC:: Phase_Zero
  PUBLIC:: Phase_Index
  PUBLIC:: Phase_Calc
  !
  TYPE:: T_Phase
  ! describes a phase, including pure phases and mixture or solution phases
    CHARACTER(LEN=23):: NamFs
    !! CHARACTER(LEN=4) :: Typ !"PURE","AQUEOUS","MIXTURE"
    !! ! "PURE" includes also "DISCRET" 
    !! ! Typ is redundant with the values of iSpc/iMix/iSol
    !! ! -> should not be included in the type declaration
    !! ! PURE     == (iSpc/=0, iMix=0,  iSol=0 )
    !! ! MIXTURE  == (iSpc=0,  iMix/=0, iSol=0 )
    !! ! SOLUTION == (iSpc=0,  iMix=0,  iSol/=0)
    INTEGER :: iSpc != index of species in vSpc, in case of pure phase 
    INTEGER :: iMix != index of mixture in vMixFas, in case of mixture phase
    INTEGER :: iSol != index of solution in vSolFas, in case of solution phase
    !
    ! variable parameters (intensive, e.g. molar, or specific)
    REAL(dp):: Grt,VolM3,WeitKg !,H,S,Cp !molar values at current (T,P), for the phase
    ! Grt=    G/RT
    ! V=      molar volume, M3/Mole
    ! WeitKg= molar weight, Kg/Mole
    !
    REAL(dp):: Mole !extensive, mole number of phase in system
    !-> should be "outside" the structure ??
  ENDTYPE T_Phase
  !
  !TYPE:: T_PhaseDat
  !  REAL(dp):: Grt,VolM3,WeitKg
  !  REAL(dp):: Mole
  !ENDTYPE T_FasData
  !
  TYPE(T_Phase):: Phase_Zero= T_Phase( &
  & "Z",          & !Name,Typ
  & 0,0,0,        & !iSpc,iMix,iSol
  & Zero,One,One, & !Grt,V,WeitKg
  & Zero)           !Mole
  !
CONTAINS

!SUBROUTINE Phase_Discrete_Init(vDisFas)
!  TYPE(T_Phase),INTENT(OUT):: vDisFas(:)
!  !
!ENDSUBROUTINE Phase_Discrete_Init

SUBROUTINE Phase_Calc(TdgK,Pbar,vSpc,vMixModel,vMixFas,Fas)
!--
!-- calculate molar properties (weight, volume, Gibbs, ...) of phase Fas
!--
  USE M_T_Species, ONLY: T_Species
  USE M_T_MixModel,ONLY: T_MixModel
  USE M_T_MixPhase,ONLY: T_MixPhase
  USE M_T_MixPhase,ONLY: MixPhase_Weit,MixPhase_GibbsRT,MixPhase_Volume
  !
  REAL(dp),        INTENT(IN)   :: TdgK,Pbar
  TYPE(T_Species), INTENT(IN)   :: vSpc(:)
  TYPE(T_MixModel),INTENT(IN)   :: vMixModel(:)
  TYPE(T_MixPhase),INTENT(IN)   :: vMixFas(:)
  TYPE(T_Phase),   INTENT(INOUT):: Fas
  !
  INTEGER :: I
  
  IF(Fas%iSpc /= 0) THEN
    !
    I=Fas%iSpc
    Fas%WeitKg= vSpc(I)%WeitKg
    Fas%Grt=    vSpc(I)%G0rt
    Fas%VolM3=  vSpc(I)%V0
    !
  ELSE IF(Fas%iMix /= 0) THEN
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
  ELSE IF (Fas%iSol /= 0) THEN
    ! TODO
  END IF
  
  !! SELECT CASE(Fas%Typ)
  !! 
  !! CASE("PURE") !,"DISCRET")
  !!   I=Fas%iSpc
  !!   
  !!   Fas%WeitKg= vSpc(I)%WeitKg
  !!   Fas%Grt=    vSpc(I)%G0rt
  !!   Fas%VolM3=  vSpc(I)%V0
  !! 
  !! CASE("MIXT") 
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
  !! END SELECT
  
  RETURN
ENDSUBROUTINE Phase_Calc

!! SUBROUTINE Phase_Zero(Fas)
!!   !
!!   TYPE(T_Phase),INTENT(OUT):: Fas
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
!! ENDSUBROUTINE Phase_Zero

INTEGER FUNCTION Phase_Index(Str,V)
!--
!-- position of phase named Str in V(1:SIZE(V))
!--
  TYPE(T_Phase),INTENT(IN)::V(:)
  CHARACTER(*), INTENT(IN)::Str
  !
  INTEGER     ::I
  
  Phase_Index=0  !IF Str not found -> I=0
  IF(SIZE(V)==0) RETURN
  
  I=0
  DO
    I=I+1
    IF(TRIM(Str)==TRIM(V(I)%NamFs)) THEN
      Phase_Index=I ; EXIT
    ENDIF
    IF(I==SIZE(V)) EXIT
  ENDDO
  
  RETURN
ENDFUNCTION Phase_Index

ENDMODULE M_T_Phase

