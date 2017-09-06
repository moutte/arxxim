module M_GEM_Vars
  !--
  !-- common vars for GEM computations
  !--
  use M_Kinds
  use M_T_Component,only: T_Component
  use M_T_MixModel, only: MaxPole
  
  implicit none

  private

  public:: GEM_Vars_Clean
  public:: vCpnGEM,tStoikioGEM
  public:: TdgK,Pbar
  public:: T_SavModel
  !
  !--types for GEM computations
  ! T_SavModel is designed to contain
  ! the different compositions that may be active
  ! for a same mixing model:
  ! For a given mixing model,
  ! there can be, in the stable assemblage,
  ! up to NPole phases of different compositions.
  ! nFas is the number of active phases for the model iModel
  type:: T_SavModel
    integer :: iModel                      ! -> mixing model
    integer :: nFas                        ! -> number of phases
    real(dp):: tXPole(1:Maxpole,1:MaxPole) ! table of phase compositions
    real(dp):: vMole(1:MaxPole)            ! mole nr of each end member
    real(dp):: vGrt0(1:MaxPole)            ! Gibbs/RT of each end member
    real(dp):: vVol0(1:MaxPole)            ! Molar Volume of each end member
  end type T_SavModel
  !
  !-- vars for GEM computations
  type(T_SavModel):: SavModelZero
  type(T_Component),allocatable:: vCpnGEM(:)
  real(dp),         allocatable:: tStoikioGEM(:,:)
  real(dp):: TdgK,Pbar

contains

subroutine GEM_Vars_Init
  SavModelZero%iModel=      0
  SavModelZero%nFas=        0
  SavModelZero%tXPole(:,:)= Zero
  SavModelZero%vMole(:)=    Zero
  SavModelZero%vGrt0(:)=    Zero
  SavModelZero%vVol0(:)=    Zero
end subroutine GEM_Vars_Init

subroutine GEM_Vars_Clean
  if(allocated(vCpnGEM))     deallocate(vCpnGEM)
  if(allocated(tStoikioGEM)) deallocate(tStoikioGEM)
end subroutine GEM_Vars_Clean

end module M_GEM_Vars

