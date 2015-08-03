MODULE M_T_SolPhase

  USE M_Kinds
  USE M_T_Species,ONLY: nElMax

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: T_SolPhase
  PUBLIC:: SolPhase_Init

  TYPE:: T_SolPhaseData
    != composition of an asymmetric solution is given
    != in terms of molalities of solute species
    REAL(dp),ALLOCATABLE:: vMolal(:)
    !
    !-- current values at (T,P,X), for the solvent and species
    REAL(dp):: LnActW                 ! solvent
    REAL(dp),ALLOCATABLE:: vLnAct(:)  ! solute species
    !
    ! REAL(dp):: LnGamW
    ! REAL(dp),ALLOCATABLE:: vLnGam(:)
    ! !
    ! !--- current values for the phase
    ! REAL(dp):: Grt,H,S,V,Cp
    !
  ENDTYPE T_SolPhaseData

  TYPE:: T_SolPhase
  != implements a phase that follows a T_SolModel mixing model,
  != i.e. an asymmetric model,
  != >> used for a molality-based solution
  !=  (e.g. electrolytes, aqueous solution)
    CHARACTER(LEN=23):: Name
    !
    !--- solution model- index in vSolModel --
    INTEGER:: iModel
    !
    !--- the composition --
    REAL(dp),ALLOCATABLE:: vXSpecies(:)
    !
  ENDTYPE T_SolPhase

CONTAINS

SUBROUTINE SolPhase_Init( & !
& iModel,    & !IN, model index
& vSolModel, & !IN, base of solution models
& S,         & !INOUT, solution phase
& Ok,Msg)      !OUT
!--
!-- initialize [solution]%iModel --
!--
  USE M_T_SolModel,ONLY: T_SolModel
  !
  INTEGER,         INTENT(IN)   :: iModel !index of solvent model
  TYPE(T_SolModel),INTENT(IN)   :: vSolModel(:)
  TYPE(T_SolPhase),INTENT(INOUT):: S
  LOGICAL,         INTENT(OUT)  :: Ok
  CHARACTER(*),    INTENT(OUT)  :: Msg
  !
  ! INTEGER:: I
  !
  Ok= .TRUE.
  Msg= ""

  ! S%iModel= 0
  ! DO I=1,SIZE(vSolModel)
  !   IF(TRIM(Str)==TRIM(vSolModel(I)%Name)) THEN
  !     S%iModel= I
  !     EXIT
  !   ENDIF
  ! ENDDO
  ! !
  ! IF(S%iModel == 0) THEN
  !   Ok= .FALSE.
  !   Msg= "Activity model "//TRIM(Str)//" Not Found"
  ! ENDIF

  S%iModel= iModel
  ALLOCATE(S%vXSpecies(vSolModel(iModel)%nSpecies))
  ! ALLOCATE(S%Dat%vMolal(vSolModel(iModel)%nSpecies))
  ! ALLOCATE(S%Dat%vLnAct(vSolModel(iModel)%nSpecies))
  !
  RETURN
END SUBROUTINE SolPhase_Init

! SUBROUTINE SolPhase_Molal_Get( & !
! & SolModel, &  !IN
! & SolFas,   &  !IN
! & nMol,     &  !IN: amount of solution phase in terms of solvent mole number!!
! & vTotCp)      !OUT
! !--
! !-- from S%vMolalCp(:), retrieve mole numbers of components /elements
! !--
!   USE M_T_SolModel,ONLY: T_SolModel
!   !
!   TYPE(T_SolModel),INTENT(IN) :: SolModel
!   TYPE(T_SolPhase),INTENT(IN) :: SolFas
!   REAL(dp),        INTENT(IN) :: nMol
!   REAL(dp),        INTENT(OUT):: vTotCp(:)
!   !
!   INTEGER:: I
!   !
!   DO I= 1,SolModel%nCp
!     vTotCp(I)= SolFas%vMolalCp(I) *nMol *SolModel%MolWeitSv
!   ENDDO
!   !
! END SUBROUTINE SolPhase_vTotCp_Get

! SUBROUTINE SolPhase_MolalCp_Set( & !
! & SolModel, &  !IN
! & vSpc,     &  !IN
! & vMol,     &  !IN
! & SolFas)      !INOUT
! !--
! !-- from a fluid compo' given in species mole nrs, set %vMolalCp(:) --
! !--
!   USE M_T_Species, ONLY: T_Species
!   USE M_T_SolModel,ONLY: T_SolModel
!   !
!   TYPE(T_SolModel),INTENT(IN)   :: SolModel
!   TYPE(T_Species), INTENT(IN)   :: vSpc(:)
!   REAL(dp),        INTENT(IN)   :: vMol(:)
!   TYPE(T_SolPhase),INTENT(INOUT):: SolFas
!   !
!   TYPE(T_Species):: S
!   INTEGER :: iSolute,iSolvent
!   INTEGER :: I,J
!   REAL(dp),ALLOCATABLE:: vX(:)
!   !
!   ALLOCATE(vX(SolModel%nCp))
!   !
!   vX(:)= Zero
!   DO I= 1,SolModel%nSpecies
!     iSolute= SolModel%vISolute(I)
!     S= vSpc(iSolute)
!     DO J=1,SolModel%nCp
!       vX(J)= vX(J) + vMol(iSolute) *S%vStoikio(J) /REAL(S%vStoikio(0))
!     ENDDO
!   ENDDO
!   !
!   iSolvent= SolModel%iSolvent
!   SolFas%vMolalCp(:)= vX(:) /vMol(iSolvent) /vSpc(iSolvent)%WeitKg
!   !
!   DEALLOCATE(vX)
!   !
! END SUBROUTINE SolPhase_MolalCp_Set

ENDMODULE M_T_SolPhase
