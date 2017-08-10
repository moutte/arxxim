module M_T_SolPhase

  use M_Kinds
  use M_T_Species,only: nElMax

  implicit none

  private

  public:: T_SolPhase
  public:: SolPhase_Init

  type:: T_SolPhaseData
    != composition of an asymmetric solution is given
    != in terms of molalities of solute species
    real(dp),allocatable:: vMolal(:)
    !
    !-- current values at (T,P,X), for the solvent and species
    real(dp):: LnActW                 ! solvent
    real(dp),allocatable:: vLnAct(:)  ! solute species
    !
    ! real(dp):: LnGamW
    ! real(dp),allocatable:: vLnGam(:)
    ! !
    ! !--- current values for the phase
    ! real(dp):: Grt,H,S,V,Cp
    !
  end type T_SolPhaseData

  type:: T_SolPhase
  != implements a phase that follows a T_SolModel mixing model,
  != i.e. an asymmetric model,
  != >> used for a molality-based solution
  !=  (e.g. electrolytes, aqueous solution)
    character(len=23):: Name
    !
    !--- solution model- index in vSolModel --
    integer:: iModel
    !
    !--- the composition --
    real(dp),allocatable:: vXSpecies(:)
    !
  end type T_SolPhase

contains

subroutine SolPhase_Init( & !
& iModel,    & !IN, model index
& nSolute,   & !IN, nSolute in vSolmodel(iModel)
& S,         & !INOUT, solution phase
& Ok,Msg)      !OUT
!--
!-- initialize [solution]%iModel --
!--
  use M_T_SolModel,only: T_SolModel
  !
  integer,         intent(in)   :: iModel !index of solvent model
  integer,         intent(in)   :: nSolute
  type(T_SolPhase),intent(inout):: S
  logical,         intent(out)  :: Ok
  character(*),    intent(out)  :: Msg
  !
  ! integer:: I
  !
  Ok= .true.
  Msg= ""

  ! S%iModel= 0
  ! do I=1,size(vSolModel)
  !   if(trim(Str)==trim(vSolModel(I)%Name)) then
  !     S%iModel= I
  !     exit
  !   end if
  ! end do
  ! !
  ! if(S%iModel == 0) then
  !   Ok= .false.
  !   Msg= "Activity model "//trim(Str)//" Not Found"
  ! end if

  S%iModel= iModel
  allocate(S%vXSpecies(nSolute))
  ! allocate(S%Dat%vMolal(vSolModel(iModel)%nSolute))
  ! allocate(S%Dat%vLnAct(vSolModel(iModel)%nSolute))
  !
  return
end subroutine SolPhase_Init

! subroutine SolPhase_Molal_Get( & !
! & SolModel, &  !IN
! & SolFas,   &  !IN
! & nMol,     &  !IN: amount of solution phase in terms of solvent mole number!!
! & vTotCp)      !OUT
! !--
! !-- from S%vMolalCp(:), retrieve mole numbers of components /elements
! !--
!   use M_T_SolModel,only: T_SolModel
!   !
!   type(T_SolModel),intent(in) :: SolModel
!   type(T_SolPhase),intent(in) :: SolFas
!   real(dp),        intent(in) :: nMol
!   real(dp),        intent(out):: vTotCp(:)
!   !
!   integer:: I
!   !
!   do I= 1,SolModel%nCp
!     vTotCp(I)= SolFas%vMolalCp(I) *nMol *SolModel%MolWeitSv
!   end do
!   !
! end subroutine SolPhase_vTotCp_Get

! subroutine SolPhase_MolalCp_Set( & !
! & SolModel, &  !IN
! & vSpc,     &  !IN
! & vMol,     &  !IN
! & SolFas)      !INOUT
! !--
! !-- from a fluid compo' given in species mole nrs, set %vMolalCp(:) --
! !--
!   use M_T_Species, only: T_Species
!   use M_T_SolModel,only: T_SolModel
!   !
!   type(T_SolModel),intent(in)   :: SolModel
!   type(T_Species), intent(in)   :: vSpc(:)
!   real(dp),        intent(in)   :: vMol(:)
!   type(T_SolPhase),intent(inout):: SolFas
!   !
!   type(T_Species):: S
!   integer :: iSolute,iSolvent
!   integer :: I,J
!   real(dp),allocatable:: vX(:)
!   !
!   allocate(vX(SolModel%nCp))
!   !
!   vX(:)= Zero
!   do I= 1,SolModel%nSolute
!     iSolute= SolModel%vISolute(I)
!     S= vSpc(iSolute)
!     do J=1,SolModel%nCp
!       vX(J)= vX(J) + vMol(iSolute) *S%vStoikio(J) /real(S%vStoikio(0))
!     end do
!   end do
!   !
!   iSolvent= SolModel%iSolvent
!   SolFas%vMolalCp(:)= vX(:) /vMol(iSolvent) /vSpc(iSolvent)%WeitKg
!   !
!   deallocate(vX)
!   !
! end subroutine SolPhase_MolalCp_Set

end module M_T_SolPhase
