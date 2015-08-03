MODULE M_GEM_Vars
  !--
  !-- common vars for GEM computations
  !--
  USE M_Kinds
  USE M_T_Component,ONLY: T_Component
  USE M_T_MixModel, ONLY: MaxPole
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: GEM_Vars_Clean
  PUBLIC:: vCpnGEM,tStoikioGEM
  PUBLIC:: TdgK,Pbar
  PUBLIC:: T_SavPhase
  !
  !--types for GEM computations
  ! T_SavPhase is designed to contain
  ! the different compositions that may be active
  ! for a same mixing model:
  ! For a given mixing model,
  ! there can be, in the stable assemblage,
  ! up to NPole phases of different compositions.
  ! nFas is the number of active phases for the model iModel
  TYPE:: T_SavPhase
    INTEGER :: iModel                      ! -> mixing model
    INTEGER :: nFas                        ! -> number of phases
    REAL(dp):: tXPole(1:Maxpole,1:MaxPole) ! table of phase compositions
    REAL(dp):: vMole(1:MaxPole)            ! mole nr of each end member
    REAL(dp):: vGrt0(1:MaxPole)            ! Gibbs/RT of each end member
    REAL(dp):: vVol0(1:MaxPole)            ! Molar Volume of each end member
  END TYPE T_SavPhase
  !
  !-- vars for GEM computations
  TYPE(T_SavPhase):: SavPhaseZero
  TYPE(T_Component),ALLOCATABLE:: vCpnGEM(:)
  REAL(dp),         ALLOCATABLE:: tStoikioGEM(:,:)
  REAL(dp):: TdgK,Pbar

CONTAINS

SUBROUTINE GEM_Vars_Init
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
END SUBROUTINE GEM_Vars_Init

SUBROUTINE GEM_Vars_Clean
  IF(ALLOCATED(vCpnGEM))     DEALLOCATE(vCpnGEM)
  IF(ALLOCATED(tStoikioGEM)) DEALLOCATE(tStoikioGEM)
ENDSUBROUTINE GEM_Vars_Clean

END MODULE M_GEM_Vars

