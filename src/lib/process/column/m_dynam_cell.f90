MODULE M_Dynam_Cell  !!_Types
  USE M_Kinds
  
  IMPLICIT NONE

  TYPE:: T_Cell_State
    !
    !! INTEGER :: nAq,nMk
    !
    REAL(dp), ALLOCATABLE:: & !
    & vMolK   (:), &          ! mole nr's of kin'phases in cell
    & vSurfK  (:), &          ! surfaces of kin'phases in cell
    & vKinQsk (:)             ! sat'states of kin'phases
    !
    CHARACTER,ALLOCATABLE:: vStatusK(:)
    !
    REAL(dp), ALLOCATABLE:: & !
    & vTotF   (:), &          !
    & vMolF   (:), &          !-> mole nr's of aqu'species in cell
    & vLnAct  (:), &          !
    & vLnGam  (:), &          !
    & vLnBuf  (:)             !
    !
    REAL(dp):: PhiF         ! vol'fraction of fluid
    !
    CONTAINS
    
    PROCEDURE :: New_
    PROCEDURE :: Free_
    PROCEDURE :: Set_
    PROCEDURE :: Get_
  
  END TYPE
  
  !! vMolK(1:nMk)=     vKinFas(1:nMk)%Dat%PhiM /vMolarVol(1:nMk) *VBox
  !! !                 !-> mole nr's of kin'phases in box
  !! vSurfK(1:nMk)=    vKinFas(1:nMk)%Dat%Surf *VBox
  !! !                 !-> surfaces of kin'phases in box
  !! vKinQsk(1:nMk)=   vKinFas(1:nMk)%Dat%QsK
  !! !                 !-> sat'states of kin'phases
  !! vStatusK(1:nMk)=  vKinFas(1:nMk)%Dat%cSat
  !! !
  !! PhiF= One - SUM(vKinFas(1:nMk)%Dat%PhiM)
  !! !                 !-> vol'fraction of fluid
  !! vMolF(1:nAq)= vSpc(vOrdAq(1:nAq))%Dat%Mole
  !! !                 !-> mole nr's of aqu'species in box
  
CONTAINS

SUBROUTINE New_( &  !
& self,          & !
& nCp,nAq,nMk    & !
& )
  !---------------------------------------------------------------------
  CLASS(T_Cell_State):: Self
  INTEGER,INTENT(IN):: nCp,nAq,nMk
  !---------------------------------------------------------------------
  ALLOCATE(               &    !
  & self%vMolK   (1:nMk), &    ! mole nr's of kin'phases in cell
  & self%vSurfK  (1:nMk), &    ! surfaces of kin'phases in cell
  & self%vKinQsk (1:nMk), &    ! sat'states of kin'phases
  & self%vStatusK(1:nMk), &    ! 
  & self%vTotF   (1:nCp), &    ! 
  & self%vMolF   (1:nAq), &    ! mole nr's of aqu'species in cell
  & self%vLnAct  (1:nAq), &    !
  & self%vLnGam  (1:nAq), &    !
  & self%vLnBuf  (1:nAq)  &    !
  & ) 

  RETURN
END SUBROUTINE New_

SUBROUTINE Set_( & !
& Self,          & !
& vMolK,         & !-> mole nr's of kin'phases in cell
& vSurfK,        & !-> surfaces of kin'phases in cell
& vKinQsk,       & !-> sat'states of kin'phases
& vStatusK,      & !
& vTotF,         & ! 
& vMolF,         & ! mole nr's of aqu'species in cell
& vLnAct,        & !
& vLnGam,        & !
& vLnBuf,        & !
& PhiF           & !-> vol'fraction of fluid
& )
  !---------------------------------------------------------------------
  CLASS(T_Cell_State):: Self
  REAL(dp)::     &    !
  & vMolK   (:), &    !-> mole nr's of kin'phases in cell
  & vSurfK  (:), &    !-> surfaces of kin'phases in cell
  & vKinQsk (:), &    !-> sat'states of kin'phases
  & vTotF   (:), &    ! 
  & vMolF   (:), &    !-> mole nr's of aqu'species in cell
  & vLnAct  (:), &    ! 
  & vLnGam  (:), &    ! 
  & vLnBuf  (:), &    ! 
  & PhiF              !-> vol'fraction of fluid
  CHARACTER:: vStatusK(:)
  !---------------------------------------------------------------------
  
  self%vMolK   (:)= vMolK   (:)
  self%vSurfK  (:)= vSurfK  (:)
  self%vKinQsk (:)= vKinQsk (:)
  self%vStatusK(:)= vStatusK(:)
  !
  self%vTotF   (:)= vTotF (:)
  self%vMolF   (:)= vMolF (:)
  self%vLnAct  (:)= vLnAct(:)
  self%vLnGam  (:)= vLnGam(:)
  !! self%vLnBuf  (:)= vLnBuf(:)
  !
  self%PhiF       = PhiF
  
  RETURN
END SUBROUTINE Set_

SUBROUTINE Free_(Self)

  CLASS(T_Cell_State):: Self

  DEALLOCATE(      &    !
  & self%vMolK   , &    ! mole nr's of kin'phases in cell
  & self%vSurfK  , &    ! surfaces of kin'phases in cell
  & self%vKinQsk , &    ! sat'states of kin'phases
  & self%vStatusK, &    ! 
  & self%vTotF,    &    ! 
  & self%vMolF,    &    ! mole nr's of aqu'species in cell
  & self%vLnAct,   &    !
  & self%vLnGam,   &    !
  & self%vLnBuf    &    !
  & ) 

  RETURN
END SUBROUTINE Free_

SUBROUTINE Get_( & !
& Self,          & !
& vMolK,         & !-> mole nr's of kin'phases in cell
& vSurfK,        & !-> surfaces of kin'phases in cell
& vKinQsk,       & !-> sat'states of kin'phases
& vStatusK,      & !
& vTotF,         & ! 
& vMolF,         & ! mole nr's of aqu'species in cell
& vLnAct,        & !
& vLnGam,        & !
& vLnBuf,        & !
& PhiF           & !-> vol'fraction of fluid
& )
  !---------------------------------------------------------------------
  CLASS(T_Cell_State):: Self
  REAL(dp)::     &    !
  & vMolK   (:), &    !-> mole nr's of kin'phases in cell
  & vSurfK  (:), &    !-> surfaces of kin'phases in cell
  & vKinQsk (:), &    !-> sat'states of kin'phases
  & vTotF   (:), &    ! 
  & vMolF   (:), &    !-> mole nr's of aqu'species in cell
  & vLnAct  (:), &    ! 
  & vLnGam  (:), &    ! 
  & vLnBuf  (:), &    ! 
  & PhiF              !-> vol'fraction of fluid
  CHARACTER:: vStatusK(:)
  !---------------------------------------------------------------------
  vMolK   (:)= self%vMolK   (:)
  vSurfK  (:)= self%vSurfK  (:)
  vKinQsk (:)= self%vKinQsk (:)
  vStatusK(:)= self%vStatusK(:)
  !
  vTotF   (:)= self%vTotF 
  vMolF   (:)= self%vMolF 
  vLnAct  (:)= self%vLnAct
  vLnGam  (:)= self%vLnGam
  !! vLnBuf  (:)= self%vLnBuf

  PhiF       = self%PhiF      
  
  RETURN
END SUBROUTINE Get_

END MODULE M_Dynam_Cell  !!_Types
