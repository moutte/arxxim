module M_Dynam_Cell  !!_Types
  use M_Kinds
  
  implicit none

  type:: T_Cell_State
    !
    !! integer :: nAq,nMk
    !
    real(dp), allocatable:: & !
    & vMolK   (:), &          ! mole nr's of kin'phases in cell
    & vSurfK  (:), &          ! surfaces of kin'phases in cell
    & vKinQsk (:)             ! sat'states of kin'phases
    !
    character,allocatable:: vStatusK(:)
    !
    real(dp), allocatable:: & !
    & vTotF   (:), &          !
    & vMolF   (:), &          !-> mole nr's of aqu'species in cell
    & vLnAct  (:), &          !
    & vLnGam  (:), &          !
    & vLnBuf  (:)             !
    !
    real(dp):: PhiF         ! vol'fraction of fluid
    !
    contains
    
    procedure :: New_
    procedure :: Free_
    procedure :: Set_
    procedure :: Get_
  
  end type
  
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
  
contains

subroutine New_( &  !
& self,          & !
& nCp,nAq,nMk    & !
& )
  !---------------------------------------------------------------------
  class(T_Cell_State):: Self
  integer,intent(in):: nCp,nAq,nMk
  !---------------------------------------------------------------------
  allocate(               &    !
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

  return
end subroutine New_

subroutine Set_( & !
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
  class(T_Cell_State):: Self
  real(dp)::     &    !
  & vMolK   (:), &    !-> mole nr's of kin'phases in cell
  & vSurfK  (:), &    !-> surfaces of kin'phases in cell
  & vKinQsk (:), &    !-> sat'states of kin'phases
  & vTotF   (:), &    ! 
  & vMolF   (:), &    !-> mole nr's of aqu'species in cell
  & vLnAct  (:), &    ! 
  & vLnGam  (:), &    ! 
  & vLnBuf  (:), &    ! 
  & PhiF              !-> vol'fraction of fluid
  character:: vStatusK(:)
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
  
  return
end subroutine Set_

subroutine Free_(Self)

  class(T_Cell_State):: Self

  deallocate(      &    !
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

  return
end subroutine Free_

subroutine Get_( & !
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
  class(T_Cell_State):: Self
  real(dp)::     &    !
  & vMolK   (:), &    !-> mole nr's of kin'phases in cell
  & vSurfK  (:), &    !-> surfaces of kin'phases in cell
  & vKinQsk (:), &    !-> sat'states of kin'phases
  & vTotF   (:), &    ! 
  & vMolF   (:), &    !-> mole nr's of aqu'species in cell
  & vLnAct  (:), &    ! 
  & vLnGam  (:), &    ! 
  & vLnBuf  (:), &    ! 
  & PhiF              !-> vol'fraction of fluid
  character:: vStatusK(:)
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
  
  return
end subroutine Get_

end module M_Dynam_Cell  !!_Types
