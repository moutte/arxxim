MODULE M_T_ChemicalSpace
!--
!-- UNUSED
!-- module for definition of a structure describing the chemical space
!--

  USE M_Kinds

  USE M_T_Element,     ONLY: T_Element
  USE M_T_Species,     ONLY: T_Species,T_SpeciesDtb
  USE M_T_MixModel,    ONLY: T_MixModel
  USE M_T_SolModel,    ONLY: T_SolModel
  USE M_T_MixPhase,    ONLY: T_MixPhase
  USE M_T_SolPhase,    ONLY: T_SolPhase
  USE M_T_Phase,       ONLY: T_Phase
  ! USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam

  USE M_T_KinFas,   ONLY: T_KinFas
  USE M_T_Kinmodel, ONLY: T_KinModel

  IMPLICIT NONE

  PRIVATE

  TYPE:: T_ChemicalSpace
    
    INTEGER:: nEle
    INTEGER:: nSpcDtb
    INTEGER:: nSpc
    INTEGER:: nMixModel
    INTEGER:: nSolModel
    INTEGER:: nMixFas
    INTEGER:: nSolFas
    INTEGER:: nFas
    ! INTEGER:: nDiscretParam
    ! INTEGER:: nDiscretModel
    
    TYPE(T_Element),    ALLOCATABLE:: vEle         (:)
    TYPE(T_SpeciesDtb), ALLOCATABLE:: vSpcDtb      (:)
    TYPE(T_Species),    ALLOCATABLE:: vSpc         (:)
    TYPE(T_MixModel),   ALLOCATABLE:: vMixModel    (:)
    TYPE(T_SolModel),   ALLOCATABLE:: vSolModel    (:)
    TYPE(T_MixPhase),   ALLOCATABLE:: vMixFas      (:)
    TYPE(T_SolPhase),   ALLOCATABLE:: vSolFas      (:)
    TYPE(T_Phase),      ALLOCATABLE:: vFas         (:)
    ! TYPE(T_DiscretParam),ALLOCATABLE:: vDiscretParam(:)
    ! TYPE(T_DiscretModel),ALLOCATABLE:: vDiscretModel(:)
      
    CONTAINS
      PROCEDURE :: New_
      PROCEDURE :: Free_
      PROCEDURE :: Set_
      PROCEDURE :: Get_
  
  END TYPE T_ChemicalSpace
  
  TYPE:: T_ChemicalState
    
    REAL(dp):: TdgK,Pbar
    
    REAL(dp),ALLOCATABLE:: vMolFracSpc(:)
    REAL(dp),ALLOCATABLE:: vMolFas(:)
  
  END TYPE T_ChemicalState
  !
  REAL(dp),ALLOCATABLE,PUBLIC:: tFormula(:,:)
  !1:nCp,1:nSp= Stoikio of All Species in Terms of Elements, in Columns
  !
  !!LOGICAL, DIMENSION(:),ALLOCATABLE,PUBLIC:: vLAqu,vLMin,vLGas 
  !a collection of masks, all of DIMENSION 1:nSp, to flag the nature of the species
  !
  TYPE(T_KinFas),   ALLOCATABLE,PUBLIC:: vKinFas(:) !-> DATA on kinetic minerals
  TYPE(T_KinModel), ALLOCATABLE,PUBLIC:: vKinModel(:) !-> DATA on kinetic models
  !
  INTEGER,PUBLIC:: &
  & nAq, &  !nAq= NrOfAqueousSpecies=  Prim'Species(1:nCp) U Second'Species(nCp+1:nAq)
  & nMn, &  !nMn= NumberOfMinerals !at the moment, nMn includes all pure phases, including gases
  & nGs !, &  !nGs= NumberOfGases !for future use 
  
  !& nCp, &  !nCp= NrOfPrimarySpecies (idem NrOfComponents) = redundant with SIZE(vCpn)
  !& nSp, &  !nSp= TotalNrOfSpecies = redundant with SIZE(vSpc)
  !& nFs, &  !number of non aqueous Multi-component phases
  !& nMk

  !TYPE:: T_GlobalState
  !  !CONTAINS a "state" of the chemical space
  !  !i.e. abundances of phases,
  !  !     abundances of species and their activities in their respective phases
  !  REAL(dp), ALLOCATABLE:: vFasMole(:)
  !  REAL(dp), ALLOCATABLE:: vSpcMole(:),vSpcLnAct(:)
  !  
  !ENDTYPE T_GlobalState
  !
CONTAINS

SUBROUTINE New_( & !
& self,          & !
& nEle, nSpcDtb,nSpc,nMixModel,nSolModel,nMixFas,nSolFas,nFas    & !
& )
  !---------------------------------------------------------------------
  CLASS(T_ChemicalSpace),INTENT(INOUT):: Self
  INTEGER:: nEle, nSpcDtb,nSpc,nMixModel,nSolModel,nMixFas,nSolFas,nFas
  !---------------------------------------------------------------------
  self%nEle      = nEle     
  self%nSpcDtb   = nSpcDtb  
  self%nSpc      = nSpc     
  self%nMixModel = nMixModel
  self%nSolModel = nSolModel
  self%nMixFas   = nMixFas  
  self%nSolFas   = nSolFas  
  self%nFas      = nFas     
  
  ! nDiscretParam
  ! nDiscretModel
  
  ALLOCATE (self%vEle         (nEle      ))
  ALLOCATE (self%vSpcDtb      (nSpcDtb   ))
  ALLOCATE (self%vSpc         (nSpc      ))
  ALLOCATE (self%vMixModel    (nMixModel ))
  ALLOCATE (self%vSolModel    (nSolModel ))
  ALLOCATE (self%vMixFas      (nMixFas   ))
  ALLOCATE (self%vSolFas      (nSolFas   ))
  ALLOCATE (self%vFas         (nFas      ))

  ! ALLOCATE (vDiscretParam(:))
  ! ALLOCATE (vDiscretModel(:))
  
  RETURN
END SUBROUTINE New_

SUBROUTINE Free_(self)
  !---------------------------------------------------------------------
  CLASS(T_ChemicalSpace),INTENT(INOUT):: Self
  !---------------------------------------------------------------------
  self%nEle      = 0
  self%nSpcDtb   = 0
  self%nSpc      = 0
  self%nMixModel = 0
  self%nSolModel = 0
  self%nMixFas   = 0
  self%nSolFas   = 0
  self%nFas      = 0
  ! nDiscretParam
  ! nDiscretModel
  !
  DEALLOCATE (self%vEle     )
  DEALLOCATE (self%vSpcDtb  )
  DEALLOCATE (self%vSpc     )
  DEALLOCATE (self%vMixModel)
  DEALLOCATE (self%vSolModel)
  DEALLOCATE (self%vMixFas  )
  DEALLOCATE (self%vSolFas  )
  DEALLOCATE (self%vFas     )
  !
  ! ALLOCATE (vDiscretParam(:))
  ! ALLOCATE (vDiscretModel(:))
  
  RETURN
END SUBROUTINE Free_

SUBROUTINE Set_( & !
& self,          & !
& vEle     ,     & !
& vSpcDtb  ,     & !
& vSpc     ,     & !
& vMixModel,     & !
& vSolModel,     & !
& vMixFas  ,     & !
& vSolFas  ,     & !
& vFas           & !
& )
  !---------------------------------------------------------------------
  CLASS(T_ChemicalSpace),INTENT(INOUT):: Self
  TYPE(T_Element),   INTENT(IN):: vEle     (:)
  TYPE(T_SpeciesDtb),INTENT(IN):: vSpcDtb  (:)
  TYPE(T_Species),   INTENT(IN):: vSpc     (:)
  TYPE(T_MixModel),  INTENT(IN):: vMixModel(:)
  TYPE(T_SolModel),  INTENT(IN):: vSolModel(:)
  TYPE(T_MixPhase),  INTENT(IN):: vMixFas  (:)
  TYPE(T_SolPhase),  INTENT(IN):: vSolFas  (:)
  TYPE(T_Phase),     INTENT(IN):: vFas     (:)
  !---------------------------------------------------------------------
  self%vEle     (:)= vEle     (:)
  self%vSpcDtb  (:)= vSpcDtb  (:)
  self%vSpc     (:)= vSpc     (:)
  self%vMixModel(:)= vMixModel(:)
  self%vSolModel(:)= vSolModel(:)
  self%vMixFas  (:)= vMixFas  (:)
  self%vSolFas  (:)= vSolFas  (:)
  self%vFas     (:)= vFas     (:)
  
  RETURN
END SUBROUTINE Set_
 
SUBROUTINE Get_( & !
& self     ,     & ! 
& vEle     ,     & !
& vSpcDtb  ,     & !
& vSpc     ,     & !
& vMixModel,     & !
& vSolModel,     & !
& vMixFas  ,     & !
& vSolFas  ,     & !
& vFas)
  !---------------------------------------------------------------------
  CLASS(T_ChemicalSpace),INTENT(IN):: Self
  TYPE(T_Element),   INTENT(OUT):: vEle     (:)
  TYPE(T_SpeciesDtb),INTENT(OUT):: vSpcDtb  (:)
  TYPE(T_Species),   INTENT(OUT):: vSpc     (:)
  TYPE(T_MixModel),  INTENT(OUT):: vMixModel(:)
  TYPE(T_SolModel),  INTENT(OUT):: vSolModel(:)
  TYPE(T_MixPhase),  INTENT(OUT):: vMixFas  (:)
  TYPE(T_SolPhase),  INTENT(OUT):: vSolFas  (:)
  TYPE(T_Phase),     INTENT(OUT):: vFas     (:)
  !---------------------------------------------------------------------
  
  vEle     (:)= self%vEle     (:)
  vSpcDtb  (:)= self%vSpcDtb  (:)
  vSpc     (:)= self%vSpc     (:)
  vMixModel(:)= self%vMixModel(:)
  vSolModel(:)= self%vSolModel(:)
  vMixFas  (:)= self%vMixFas  (:)
  vSolFas  (:)= self%vSolFas  (:)
  vFas     (:)= self%vFas     (:)
  
  RETURN
END SUBROUTINE Get_

ENDMODULE M_T_ChemicalSpace


