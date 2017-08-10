module M_Global_Alloc
  use M_Kinds
  use M_Trace,  only: fTrc,T_,Stop_,iDebug
  implicit none
  !
  private
  !
  public:: Elements_Alloc_forDtb
  public:: SpeciesDtb_Alloc
  public:: MixModels_Alloc
  public:: MixPhases_Alloc
  public:: Phases_Alloc
  public:: Phases_Alloc_New
  public:: DiscretParam_Alloc
  !
contains

subroutine Elements_Alloc_forDtb(vEle,nEle)
!--
!-- read all elements (name,valency, weight) from database -> builds vEle
!--
  use M_T_Element,   only: T_Element,Element_Index
  use M_Element_Read,only: T_LnkEle,Elements_BuildLnk,Elements_LnkToVec
  !
  type(T_Element),allocatable,intent(out):: vEle(:)
  integer,                    intent(out):: nEle
  !
  type(T_LnkEle),pointer:: LnkEle
  !
  if(idebug>1) write(fTrc,"(/,A)") "< Elements_Alloc_forDtb"
  !
  call Elements_BuildLnk(LnkEle,nEle)
  !
  if(nEle>0) then
    !
    if(allocated(vEle)) deallocate(vEle)
    allocate(vEle(1:nEle))
    !
    call Elements_LnkToVec(LnkEle,vEle)
    !
    if(Element_Index("O__",vEle)==0) call Stop_("Oxygen not Found ???") !________stop 
    if(Element_Index("H__",vEle)==0) call Stop_("Hydrogen not Found ???") !______stop 
    !
  end if
  !
  if(idebug>1) write(fTrc,"(A,/)") "</ Elements_Alloc_forDtb"
end subroutine Elements_Alloc_forDtb

subroutine SpeciesDtb_Alloc(vSpcDtb,N) !-> build vSpcDtb
  use M_T_Species,      only: T_SpeciesDtb
  use M_SpeciesDtb_Read,only: T_LnkSpc,SpeciesDtb_LnkToVec,SpeciesDtb_BuildLnk
  !
  type(T_SpeciesDtb),allocatable,intent(out):: vSpcDtb(:)
  integer,                       intent(out):: N
  !
  type(T_LnkSpc), pointer:: LnkSpc
  !
  call SpeciesDtb_BuildLnk(LnkSpc,N)
  !-> must come after Dtb_Read_HSV, which builds all vDtbXxx from HSV databases
  !-> from vDtbXxxXxx tables, build linked list LnkSpc
  !
  if(N>0) then !transfer LnkSpc to vSpcDtb
    if(allocated(vSpcDtb)) deallocate(vSpcDtb)
    allocate(vSpcDtb(1:N))
    call SpeciesDtb_LnkToVec(LnkSpc,vSpcDtb)
  end if
  !
  !! if(iDebug>3) then 
  !!   allocate(tFormul(1:size(vEle),1:size(vSpc)))
  !!   call FormulaTable_Calc(vEle,vSpc,tFormul)
  !!   call FormulaTable_Sho(vEle,vSpc,tFormul)
  !!   deallocate(tFormul)
  !! end if
  !
end subroutine SpeciesDtb_Alloc

subroutine MixModels_Alloc(vSpc,    vMixModel) !-> build vMixModel
  use M_T_Species,    only: T_Species
  use M_T_MixModel,   only: T_MixModel
  use M_MixModel_Read,only: T_LnkMixModel,MixModel_BuildLnk,MixModel_LnkToVec
  !
  type(T_Species), dimension(:),            intent(in) :: vSpc
  type(T_MixModel),dimension(:),allocatable,intent(out):: vMixModel
  !
  type(T_LnkMixModel),pointer:: LnkSol
  integer::N
  logical:: Ok
  character(len=80):: MsgError
  !
  call MixModel_BuildLnk(vSpc,LnkSol,N,Ok,MsgError)
  !
  if(.not.Ok) then
    call Stop_("Error reading MIXTURE.MODEL: "//trim(MsgError))
  end if
  !
  if(N>0) then
    if(allocated(vMixModel)) deallocate(vMixModel)
    allocate(vMixModel(1:N))
    call MixModel_LnkToVec(LnkSol,vMixModel)
  end if
  !
end subroutine MixModels_Alloc

subroutine MixPhases_Alloc(vSpc,vMixModel,    vMixFas)
  use M_T_Species,    only: T_Species
  use M_T_MixModel,   only: T_MixModel
  use M_T_MixPhase,   only: T_MixPhase
  use M_T_MixPhase,   only: MixPhase_Zero
  use M_MixPhase_Read,only: T_LnkFas,MixPhase_BuildLnk,MixPhase_LnkToVec
  !
  type(T_Species),             intent(in)   :: vSpc(:)
  type(T_MixModel),            intent(in)   :: vMixModel(:)
  type(T_MixPhase),allocatable,intent(inout):: vMixFas(:)
  !
  type(T_LnkFas),pointer:: Lnk
  integer::I,nMixFas
  logical:: Ok
  character(len=80):: MsgError
  !
  if(idebug>1) write(fTrc,'(/,A)') "< MixPhases_Alloc"
  !
  call MixPhase_BuildLnk(vSpc,vMixModel,nMixFas,Lnk,Ok,MsgError)
  !
  if(.not.Ok) then
    call Stop_("Error reading MIXTURE block, "//trim(MsgError))
  end if
  !
  if(nMixFas>0) then
    !
    if(allocated(vMixFas)) deallocate(vMixFas)
    allocate(vMixFas(1:nMixFas))
    !
    do I=1,nMixFas; call MixPhase_Zero(vMixFas(I)); end do
    !
    call MixPhase_LnkToVec(Lnk,vMixFas)
    !
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ MixPhases_Alloc"
end subroutine MixPhases_Alloc

subroutine Phases_Alloc(vSpc,vMixFas,    vFas)
!--
!-- -> build vFas
!-- vFas- (H2O) U (non-aqueous pure species) U (solution phases)
!--
  use M_T_Phase,    only: T_Phase
  use M_T_Species,  only: T_Species,Species_Index
  use M_T_MixPhase, only: T_MixPhase
  use M_T_Phase,    only: Phase_Zero
  !
  type(T_Species),          intent(in) :: vSpc(:)
  type(T_MixPhase),         intent(in) :: vMixFas(:)
  type(T_Phase),allocatable,intent(out):: vFas(:)
  !
  integer:: I,N
  !
  ! calculate number of non-aqueous species -> included as "PURE" phase
  N=0
  if(Species_Index("H2O",vSpc)/=0) N=N+1
  ! "H2O" is "AQU" but it should be included here as pure
  N= N &
  & + count(vSpc(:)%Typ/="AQU") !& !skip the solute species
  !& + size(vMixFas)
  !
  if(allocated(vFas)) deallocate(vFas)
  allocate(vFas(1:N)); vFas(:)=Phase_Zero
  !
  N=0
  !
  if(Species_Index("H2O",vSpc)/=0) then
    ! include H2O species as pure phase in phase table
    N=N+1
    vFas(N)%NamFs= "H2O"
    !! vFas(N)%Typ=   "PURE"
    vFas(N)%iSpc=  Species_Index("H2O",vSpc)
    vFas(N)%iMix=  0
    vFas(N)%iSol=  0
  end if
  !
  do I=1,size(vSpc)
    ! include all non'aqu thermodyn'species and discretized species
    ! as pure phases in phase table
    if(vSpc(I)%Typ/="AQU") then
      N=N+1
      vFas(N)%NamFs= trim(vSpc(I)%NamSp)
      !! vFas(N)%Typ=  "PURE"
      !! if(vSpc(I)%iDiscret>0) vFas(N)%Typ= "DISCRET"
      vFas(N)%iSpc= I
      vFas(N)%iSol= 0
      vFas(N)%iMix= 0
    end if
  end do
  !
  ! do I=1,size(vMixFas)
  !   ! include all non'aqu solution phases in phase table
  !   vFas(N+I)%NamFs= trim(vMixFas(I)%Name)
  !   vFas(N+I)%Typ=   "MIXT"
  !   vFas(N+I)%iSol=  0
  !   vFas(N+I)%iSpc=  0
  !   vFas(N+I)%iMix=  I
  ! end do
  
  ! write(11,'(A)') "<==Phases_Alloc====="
  ! do i=1,size(vFas)
  !   write(11,'(I3,1X,A)') i, vFas(i)%NamFs
  ! end do
  ! write(11,'(A)') "===================>"
  
  return
end subroutine Phases_Alloc

subroutine Phases_Alloc_New(vSpc,vMixFas,vSolFas,    vFas)
!--
!-- -> build vFas
!-- vFas- (non-aqueous pure species) U (mixture phases) U (solution phases)
!--
  use M_T_Species,  only: T_Species,Species_Index
  use M_T_SolPhase, only: T_SolPhase
  use M_T_MixPhase, only: T_MixPhase
  use M_T_Phase,    only: Phase_Zero,T_Phase
  !
  type(T_Species),          intent(in) :: vSpc(:)
  type(T_MixPhase),         intent(in) :: vMixFas(:)
  type(T_SolPhase),         intent(in) :: vSolFas(:)
  type(T_Phase),allocatable,intent(out):: vFas(:)
  !
  integer:: I,N
  !
  ! calculate number of non-aqueous species -> included as "PURE" phase
  N= count(vSpc(:)%Typ/="AQU") &
  &  + size(vMixFas) &
  &  + size(vSolFas)
  !
  !! ! calculate number of non-aqueous species -> included as "PURE" phase
  !! N=0
  !! if(Species_Index("H2O",vSpc)/=0) N=N+1
  !! ! "H2O" is "AQU" but it should be included here as pure
  !! N= N &
  !! & + count(vSpc(:)%Typ/="AQU") & !skip the solute species
  !! & + size(vMixFas)
  !
  if(allocated(vFas)) deallocate(vFas)
  allocate(vFas(1:N)); vFas(:)=Phase_Zero
  !
  N=0
  !
  !! if(Species_Index("H2O",vSpc)/=0) then
  !! ! include H2O species as pure phase in phase table
  !!   N=N+1
  !!   vFas(N)%NamFs= "H2O"
  !!   vFas(N)%Typ=   "PURE"
  !!   vFas(N)%iSpc=  Species_Index("H2O",vSpc)
  !! end if
  !
  do I=1,size(vSpc)
    ! include all non'aqu thermodyn'species (and discretized species)
    ! as pure phases in phase table
    if(vSpc(I)%Typ/="AQU") then
      N= N+1
      vFas(N)%NamFs= trim(vSpc(I)%NamSp)
      !! vFas(N)%Typ=  "PURE"
      ! nota:
      !   discretized phases are comprised within the "pure phases"
      vFas(N)%iSpc= I
      vFas(N)%iSol= 0
      vFas(N)%iMix= 0
    end if
  end do
  !
  do I=1,size(vSolFas)
    ! include all non'aqu solution phases in phase table
    N= N+1
    vFas(N)%NamFs= trim(vSolFas(I)%Name)
    !! vFas(N)%Typ=  "SOLU"
    vFas(N)%iSol= I
    vFas(N)%iSpc= 0
    vFas(N)%iMix= 0
  end do
  !
  do I=1,size(vMixFas)
    ! include all non'aqu solution phases in phase table
    N= N+1
    vFas(N)%NamFs= trim(vMixFas(I)%Name)
    !! vFas(N)%Typ=  "MIXT"
    vFas(N)%iSol= 0
    vFas(N)%iSpc= 0
    vFas(N)%iMix= I
  end do
  
  ! write(11,'(A)') "<==Phases_Alloc_New=="
  ! do i=1,size(vFas)
  !   write(11,'(I3,1X,A)') i, vFas(i)%NamFs
  ! end do
  ! write(11,'(A)') "====================>"
  
  return
end subroutine Phases_Alloc_New

subroutine DiscretParam_Alloc(vDiscretModel, vDiscretParam)
!--
!-- build vDiscretParam, array of T_DiscretParam
!-- that will be used to build pseudo-species derived from discretized mixtures(s)
!--
  use M_T_MixModel,only: T_MixModel
  use M_T_DiscretModel !T_DiscretModel,T_DiscretParam
  !
  type(T_DiscretModel),            intent(in)   :: vDiscretModel(:)
  type(T_DiscretParam),allocatable,intent(inout):: vDiscretParam(:)
  !
  integer:: N,iM
  !
  if(idebug>1) write(fTrc,'(/,A)') "< DiscretParam_Alloc"
  !
  N= 0
  do iM= 1, size(vDiscretModel)
    N= N + vDiscretModel(iM)%DimTot
  end do
  deallocate(vDiscretParam)
  allocate(vDiscretParam(N))
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ DiscretParam_Alloc"
  !
end subroutine DiscretParam_Alloc

end module M_Global_Alloc

