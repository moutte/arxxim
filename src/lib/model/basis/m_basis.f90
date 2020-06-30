 module M_Basis
!--
!-- basic tools used by higher level modules:
!-- build vFas, vMixModel, ...
!-- operations on basis
!--
  use M_Kinds
  use M_Trace,only: fTrc,T_,Stop_,iDebug
  !
  implicit none
  !
  private
  !
  public:: Basis_Build
  public:: Basis_Change
  public:: Basis_Change_Wrk
  public:: Basis_InitBalance
  public:: Basis_CpnInert
  !
contains

subroutine Basis_Build(LBuffer,  vCpn) !,nCi,nAx,nMx,nAs,nMs)
  use M_T_Component,only: T_Component
  use M_Basis_Tools
  !
  use M_Global_Vars,only: vEle,vSpc,vSpcDat,vMixModel
  use M_Global_Vars,only: vFas,vMixFas
  use M_Global_Vars,only: tFormula
  !
  use M_System_Vars,only: iO_,iH_,iOx,BufferIsExtern
  !
  use M_Basis_Vars,only: iBal
  use M_Basis_Vars,only: tStoikio
  use M_Basis_Vars,only: vLCi,vLAx,vLMx,vLAs,vLMs
  use M_Basis_Vars,only: nCi, nAx, nMx, nAs, nMs
  use M_Basis_Vars,only: vOrdPr,vOrdAq,vOrdAs,vOrdMs
  use M_Basis_Vars,only: vPrmFw,vPrmBk
  use M_Basis_Vars,only: tAlfSp,tAlfPr,tAlfAs,tAlfMs
  use M_Basis_Vars,only: tNuSp,tNuAs,tNuMs
  use M_Basis_Vars,only: tNuFas,tAlfFs
  !---------------------------------------------------------------------
  logical,          intent(in)   :: LBuffer
  type(T_Component),intent(inout):: vCpn(:)
  !integer,          intent(out)  :: nCi, nAx, nMx, nAs, nMs 
  !---------------------------------------------------------------------
  integer:: nSp,nCp,nFs,I
  integer:: nAq,nMn
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Basis_Build"
  !
  nSp= size(vSpc)
  nCp= size(vEle)
  nFs= size(vFas)
  !
  nAq= count(vSpc%Typ(1:3)=="AQU")
  nMn= count(vSpc%Typ(1:3)=="MIN")+count(vSpc%Typ(1:3)=="GAS")
  !
  call Basis_Dimensions(   & !
  & BufferIsExtern,        & !IN
  & vCpn,vSpc,nAq,nMn,     & !IN
  & nCi,nMx,nAx,nMs,nAs)     !OUT, may change when components'status change
  !
  call Basis_Alloc(nSp,nCp,nFs,nAq,nAs,nMs)
  !
  !------------------------------------------- build logical tables vLxx
  call Basis_SpcTypes_Build( &
  & vSpc,vCpn, &                 !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs)    !OUT
  !
  !-------------------------------------- order components, inert/mobile
  call Basis_Cpn_Order( &
  & vSpc,count(vLCi),count(vLAx), & !IN
  & vSpcDat,vCpn)                   !INOUT
  !
  !-------------------------------- indexes of elements O, H, Ox in vCpn
  !------ (may change in basis change, because element'order may change)
  call Basis_Calc_iO_iH_iOx( &
  & vEle,vCpn, &
  & iO_,iH_,iOx)
  !
  !------------ component status and order may have changed -> find iBal
  if(all(vCpn(:)%Statut=="INERT")) then
    iBal=0 !-> no balance element -> mass balance on H
  else
    call Basis_FindIBal(vCpn,iH_,    iBal)
    if(iBal==0) iBal= iH_
    ! when iBal=iH_, mass balance on H is "replaced" by charge balance
  end if
  !
  !--------------- component order may have changed -> new stoikio table
  tStoikio(1:nCp,1:nSp)= tFormula(vCpn(1:nCp)%iEle,1:nSp) !-> line permutation
  !
  !------------------------------------------------- build vOrdXx,vPrmXx
  call Basis_SpcOrdPrm_Build(    & !
  & vSpc,vCpn,                   & !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs,    & !IN
  & vOrdPr,vOrdAq,vOrdAs,vOrdMs, & !OUT
  & vPrmFw,vPrmBk)                 !OUT
  !
  if(iDebug>1) &
  & call Basis_StoikTable_Sho( &
  & vEle,vCpn,vSpc,tStoikio,   & !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs)    !IN                                                                          
  !
  !-------------------------------------------------- build tAlfXx,tNuXx
  call Basis_AlfaNu_Build( &
  & vCpn,vSpc,tStoikio,          & !IN
  & vFas,vMixFas,vMixModel,      & !
  & vLCi,vLAx,vLMx,vLAs,vLMs,    & !
  & iDebug>1,LBuffer,            & !
  & vOrdPr,vOrdAs,vOrdMs,        & !
  & vCpn(:)%Mole,                & !INOUT
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !OUT
  & tNuSp,tNuAs,tNuMs,           & !
  & tNuFas,tAlfFs)                 !
  !
  !------------------------------------- modify stoichio tables for iBal
  if(iBal/=0) call Basis_InitBalance( & !
  & iBal,                 & !IN
  & vOrdPr,vOrdAs,vOrdMs, & !IN
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !INOUT
  & tAlfFs)                 !INOUT
  !
  do I=1,nCp
    if(vCpn(I)%Statut/="INERT") vCpn(iBal)%Mole= Zero
  end do
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Basis_Build"
  !
end subroutine Basis_Build

subroutine Basis_Change(Cod,    vCpn) !,nCi,nAx,nMx,nAs,nMs)
  !
  use M_T_Component,only: T_Component
  ! use M_T_Species,  only: T_SpcData
  use M_Dtb_Const,  only: T_CK
  use M_Basis_Tools,only: Basis_Cpn_ISpc_Select,Basis_MixPhase_CpnUpdate
  !
  use M_Global_Vars,only: vSpc,vSpcDat,vMixFas,SolModel
  use M_Basis_Vars, only: tStoikio
  !---------------------------------------------------------------------
  character(len=3), intent(in)   :: Cod
  ! type(T_SpcData),  intent(in)   :: vSpcDat(:)
  type(T_Component),intent(inout):: vCpn(:)
  !integer,          intent(out)  :: nCi,nAx,nMx,nAs,nMs
  !---------------------------------------------------------------------
  integer:: nCp,I,nAq
  logical:: LBuffer
  integer,allocatable:: vIndex(:)
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Basis_Change"
  !
  LBuffer= (Cod(1:3)=="DYN") !! .true. !! !.or. (Cod(1:2)=="EQ")
  !
  nCp= size(vCpn)
  nAq= count(vSpc%Typ(1:3)=="AQU")
  !
  if(size(vMixFas)>0) &
  & call Basis_MixPhase_CpnUpdate( &
  & vMixFas,   &  !in
  & vSpc,      &  !in
  & vSpcDat,vCpn) !inout
  !
  select case(trim(Cod))
  !
  case("DYN","EQU")
    !
    where(vCpn(:)%Statut/="BUFFER") vCpn(:)%Statut="INERT"
    !
    ! print *,"sum(tStoikio(I,:) *vSpcDat(:)%Mole)"
    ! do I=1,nCp
    !   vCpn(I)%Mole= sum(tStoikio(I,:) *vSpcDat(:)%Mole)
    !   print '(I3,G15.6,A)',I,vCpn(I)%Mole,trim(vEle(vCpn(I)%iEle)%NamEl)
    ! end do
    !
    ! compute element mole numbers refer to abundances in fluid only
    ! (-> normal case when coming from a speciation calc')
    !print *,"sum(tStoikio(I,:) *vSpcDat(:)%Mole)"
    do I=1,size(vCpn)
      vCpn(I)%Mole= &
      & sum(tStoikio(I,1:nAq)*vSpcDat(1:nAq)%Mole)
      !-> abundance in fluid
      !print '(I3,G15.6,A)',I,vCpn(I)%Mole,trim(vEle(vCpn(I)%iEle)%NamEl)
      ! & + sum(tAlfFs(I,1:nFs) *vFas(1:nFs)%Mole)
      ! write(FF,'(A,A1,G15.6)') vEle(vCpn(I)%iEle)%NamEl,T_,X
    end do
    ! call pause_
    !
    vCpn(1:nCp)%LnAct= vSpcDat(vCpn(1:nCp)%iSpc)%LAct
    !
    !----------------------- assign prim'species to each inert component
    allocate(vIndex(nCp))
    !
    call Basis_Cpn_ISpc_Select( & !
    & SolModel,vSpc,vSpcDat,vCpn, & !in
    & vIndex) !inout, only values of vCpn(:)%iSpc are possibly changed
    !
    vCpn(:)%iSpc= vIndex(:)
    deallocate(vIndex)
    !----------------------/ assign prim'species to each inert component
    !
    !! if(iOx/=0) call Basis_RedoxAssign(vCpn,vEle,vSpc)
  end select
  !
  call Basis_Build(LBuffer,vCpn) !,nCi,nAx,nMx,nAs,nMs)
  !-> include following
  ! Basis_Dimensions      !-> nCi,nMx,nAx,nMs,nAs
  ! Basis_Alloc           !-> allocate basis_vars
  ! Basis_SpcTypes_Build  !-> build vLCi,vLAx,vLMx,vLAs,vLMs
  ! Basis_Cpn_Order       !-> order components, inert/mobile

  ! Basis_Calc_iO_iH_iOx  !-> indexes of elements O, H, Ox in vCpn
  ! Basis_SpcOrdPrm_Build !-> vOrdPr,vOrdAq,vOrdAs,vOrdMs,vPrmFw,vPrmBk
  ! Basis_AlfaNu_Build    !-> build tAlfXx,tNuXx
  ! Basis_InitBalance     !
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Basis_Change"
  !
end subroutine Basis_Change

subroutine Basis_Change_Wrk(vCpn) !,nCi,nAx,nMx,nAs,nMs)
  !---------------------------------------------------------------------
  use M_T_Component,only: T_Component
  use M_Basis_Tools,only: Basis_Cpn_ISpc_Select
  !
  use M_Global_Vars,only: vSpc,vSpcDat,SolModel
  ! 
  type(T_Component),intent(inout):: vCpn(:)
  !integer,          intent(out)  :: nCi,nAx,nMx,nAs,nMs
  !---------------------------------------------------------------------
  integer,allocatable:: vIndex(:)
  !---------------------------------------------------------------------
  if(iDebug==4) write(fTrc,'(/,A)') "< Basis_Change_wrk"
  !
  allocate(vIndex(size(vCpn)))
  !
  call Basis_Cpn_ISpc_Select( &
  & SolModel,vSpc,vSpcDat,vCpn, & !in
  & vIndex)
  !
  !------------ change basis only if the primary species set has changed
  if(any ( vCpn(:)%iSpc /= vIndex(:)) ) then
    !
    vCpn(:)%iSpc= vIndex(:)
    !
    call Basis_Build(.false.,vCpn)
    !
    if(iDebug>2) write(*,'(A)') "---------------< basis changed >----"
    !
  end if
  !--//
  deallocate(vIndex)
  !
  if(iDebug==4) write(fTrc,'(A,/)') "</ Basis_Change_wrk"
end subroutine Basis_Change_Wrk

subroutine Basis_IndexInit(vSpc,nAq,nMn) ! (isH_,isOH,isO2)
!--
!-- calc' dimensions and indexes that are constant for given vSpc
!--
  use M_T_Species,only: T_Species
  !---------------------------------------------------------------------
  type(T_Species), intent(in) :: vSpc(:)
  integer,         intent(out):: nAq,nMn
  !---------------------------------------------------------------------
  ! integer, intent(out):: isH_,isOH,isO2
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Basis_IndexInit"
  !
  nAq= count(vSpc%Typ(1:3)=="AQU")
  nMn= count(vSpc%Typ(1:3)=="MIN") + count(vSpc%Typ(1:3)=="GAS")
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Basis_IndexInit"
  !
end subroutine Basis_IndexInit

subroutine Basis_InitBalance( & !
& iBal,                 & !IN
& vOrdPr,vOrdAs,vOrdMs, & !IN
& tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs) !INOUT
!--
!-- replace current stoichio value by the species charge
!-- in line iBal of all stoikio matrices tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs
!--
  !---------------------------------------------------------------------
  use M_Global_Vars,only: vSpc
  !---------------------------------------------------------------------
  integer, intent(in)   :: iBal
  integer, intent(in)   :: vOrdPr(:),vOrdAs(:),vOrdMs(:)
  real(dp),intent(inout):: tAlfSp(:,:),tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:),tAlfFs(:,:)
  !---------------------------------------------------------------------
  tAlfSp(iBal,:)=                    real(vSpc(:)%Z)
  tAlfPr(iBal,:)=                    real(vSpc(vOrdPr(:))%Z)
  if(size(vOrdAs)>0) tAlfAs(iBal,:)= real(vSpc(vOrdAs(:))%Z)
  if(size(vOrdMs)>0) tAlfMs(iBal,:)= real(vSpc(vOrdMs(:))%Z)
  tAlfFs(iBal,:)= Zero
  !
end subroutine Basis_InitBalance

subroutine Basis_CpnInert(tStoikio,vCpnIn,vSpcDat,vCpnInert)
!--
!-- used only in Equil_Write,
!-- for printing composition of aqueous phase as closed system
!-- make "INERT" all "MOBILE" components,
!-- reorder components, assign prim'species, ...
!--
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_SpcData
  ! use M_Basis_Tools,only: Basis_Species_Select_Test
  use M_Basis_Tools,only: Basis_Cpn_ISpc_Select
  !
  use M_Global_Vars,only: vSpc,SolModel
  !
  real(dp),         intent(in) :: tStoikio(:,:)
  type(T_Component),intent(in) :: vCpnIn(:)
  type(T_SpcData),  intent(in) :: vSpcDat(:)
  type(T_Component),intent(out):: vCpnInert(:)
  !
  integer:: i
  integer,allocatable:: vIndex(:)
  !
  vCpnInert= vCpnIn
  vCpnInert(:)%Statut="INERT"
  !
  !----------- compute total amount of each component in the fluid phase
  do i=1,size(vCpnInert)
    vCpnInert(i)%Mole= &
    & sum( tStoikio(i,:)*vSpcDat(:)%Mole, mask= (vSpc(:)%Typ=="AQU") )
  end do
  !
  !------------------------------------------------- assign prim'species
  allocate(vIndex(size(vCpnInert)))
  !
  call Basis_Cpn_ISpc_Select( & !
  & SolModel,vSpc,vSpcDat,vCpnInert, & !in
  & vIndex) !out
  !
  vCpnInert(:)%iSpc= vIndex(:)
  !------------------------------------------------/ assign prim'species
  !
  deallocate(vIndex)
  !
end subroutine Basis_CpnInert

end module M_Basis


