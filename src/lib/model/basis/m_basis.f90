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

subroutine Basis_Build(LBuffer,vCpn)
  use M_T_Component,only: T_Component
  use M_Basis_Tools
  !
  use M_Global_Vars,only: vEle,vSpc,vFas,vMixFas,vMixModel
  use M_Global_Vars,only: tFormula
  use M_Global_Vars,only: nAq,nMn,nGs
  use M_System_Vars,only: BufferIsExtern
  !
  use M_Basis_Vars,only: iO_,iH_,iOx,iBal
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
  !---------------------------------------------------------------------
  integer:: nSp,nCp,nFs,I
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Build"
  !
  nSp= size(vSpc)
  nCp= size(vEle)
  nFs= size(vFas)
  !
  call Basis_Dimensions(   & !
  & BufferIsExtern,        & !IN
  & vCpn,vSpc,nAq,nMn,nGs, & !IN
  & nCi,nMx,nAx,nMs,nAs)     !OUT, may change when components'status change
  !
  call Basis_Alloc(nSp,nCp,nFs,nAq,nAs,nMs)
  !
  !---------------------------------------- build logical tables vLxx --
  call Basis_SpcTypes_Build( &
  & vSpc,vCpn, &                 !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs)    !OUT
  !
  !----------------------------------- order components, inert/mobile --
  call Basis_Cpn_Order( &
  & vSpc,count(vLCi),count(vLAx), & !IN
  & vCpn) !INOUT
  !
  !----------------------------- indexes of elements O, H, Ox in vCpn --
  !--- (may change in basis change, because element'order may change) --
  call Basis_Calc_iO_iH_iOx( &
  & vEle,vCpn, &
  & iO_,iH_,iOx)
  !
  !--------- component status and order may have changed -> find iBal --
  if(all(vCpn(:)%Statut=="INERT")) then
    iBal=0 !-> no balance element -> mass balance on H
  else
    call Basis_FindIBal(vCpn,vSpc,iH_,iBal)
    if(iBal==0) iBal= iH_
    ! when iBal=iH_, mass balance on H is "replaced" by charge balance
  end if
  !
  !------------ component order may have changed -> new stoikio table --
  tStoikio(1:nCp,1:nSp)= tFormula(vCpn(1:nCp)%iEle,1:nSp) !-> line permutation
  !
  ! !--- NEW2010-10-29
  ! if(CpnIsSpc) then
  !   do i=1,nCp
  !     !vCpn(i)%iEle= I
  !     vCpn(i)%vStoikCp(0:nCp+1)= vSpc(vCpn(i)%iSpc)%vStoikio(0:nCp+1)
  !   end do
  ! else
  !   do i=1,nCp
  !     vCpn(i)%iEle= I
  !     vCpn(i)%vStoikCp(0)= 1 !formula divider !!
  !     vCpn(i)%vStoikCp(:)= 0
  !     vCpn(i)%vStoikCp(i)= 1
  !   end do
  ! end if
  ! !---/ NEW2010-10-29
  !
  !---------------------------------------------- build vOrdXx,vPrmXx --
  call Basis_SpcOrdPrm_Build(    & !
  & vSpc,vCpn,                   & !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs,    & !IN
  & vOrdPr,vOrdAq,vOrdAs,vOrdMs, & !OUT
  & vPrmFw,vPrmBk)                 !OUT
  !
  if(iDebug>1) &
  & call Basis_StoikTable_Sho(vEle,vCpn,vSpc,tStoikio,vLCi,vLAx,vLMx,vLAs,vLMs)
  !
  !----------------------------------------------- build tAlfXx,tNuXx --
  call Basis_AlfaNu_Build( &
  & vEle,vCpn,vSpc,tStoikio,     & !IN
  & vFas,vMixFas,vMixModel,      & !
  & vLCi,vLAx,vLMx,vLAs,vLMs,    & !
  & iDebug>1,LBuffer,            & !
  & vOrdPr,vOrdAs,vOrdMs,        & !
  & vCpn(:)%Mole,                & !INOUT
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !OUT
  & tNuSp,tNuAs,tNuMs,           & !
  & tNuFas,tAlfFs)                 !
  !
  !---------------------------------- modify stoichio tables for iBal --
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
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Build"
  !
end subroutine Basis_Build

subroutine Basis_Change( & !
& Cod,                   & !
!! & tStoikio,           & !
!! & isW,isH_,isOH,isO2, & !
!! & MWsv,               & !
& vCpn)
  
  use M_T_Component,only: T_Component
  use M_Dtb_Const,  only: T_CK
  use M_Basis_Tools
  !
  use M_Global_Vars,only: vSpc,nAq,vMixFas
  use M_Basis_Vars, only: isW,isH_,isOH,isO2,tStoikio,MWsv
  !
  character(len=3), intent(in)   :: Cod
  !! real(dp),         intent(in)   :: tStoikio(:,:)
  !! integer,          intent(out)  :: isW,isH_,isOH,isO2
  !! real(dp),         intent(out)  :: MWsv
  type(T_Component),intent(inout):: vCpn(:)
  !
  integer:: nCp,I
  logical:: LBuffer
  integer,allocatable:: vIndex(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Change"
  !
  call Basis_IndexInit(isW,isH_,isOH,isO2,MWsv)
  !
  LBuffer= (Cod(1:3)=="DYN") !! .true. !! !.or. (Cod(1:2)=="EQ")
  !
  nCp= size(vCpn)
  !
  if(size(vMixFas)>0) &
  & call Basis_MixPhase_CpnUpdate( &
  & vMixFas,   & !in
  & vSpc,vCpn)   !inout
  !
  select case(trim(Cod))
  !
  case("DYN","EQU")
    !
    where(vCpn(:)%Statut/="BUFFER") vCpn(:)%Statut="INERT"
    !
    ! print *,"SUM(tStoikio(I,:) *vSpc(:)%Dat%Mole)"
    ! do I=1,nCp
    !   vCpn(I)%Mole= SUM(tStoikio(I,:) *vSpc(:)%Dat%Mole)
    !   print '(I3,G15.6,A)',I,vCpn(I)%Mole,trim(vEle(vCpn(I)%iEle)%NamEl)
    ! end do
    !
    ! compute element mole numbers refer to abundances in fluid only
    ! (-> normal case when coming from a speciation calc')
    !print *,"SUM(tStoikio(I,:) *vSpc(:)%Dat%Mole)"
    do I=1,size(vCpn)
      vCpn(I)%Mole= &
      & SUM(tStoikio(I,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> abundance in fluid
      !print '(I3,G15.6,A)',I,vCpn(I)%Mole,trim(vEle(vCpn(I)%iEle)%NamEl)
      ! & + SUM(tAlfFs(I,1:nFs) *vFas(1:nFs)%Mole)
      ! write(FF,'(A,A1,G15.6)') vEle(vCpn(I)%iEle)%NamEl,T_,X
    end do
    ! call pause_
    !
    vCpn(1:nCp)%LnAct= vSpc(vCpn(1:nCp)%iSpc)%Dat%LAct
    !
    !-------------------- assign prim'species to each inert component --
    allocate(vIndex(nCp))
    !
    call Basis_Cpn_ISpc_Select( & !
    & vSpc,isW,vCpn, & !in
    & vIndex) !inout, only values of vCpn(:)%iSpc are possibly changed
    !
    vCpn(:)%iSpc= vIndex(:)
    deallocate(vIndex)
    !-------------------/ assign prim'species to each inert component --
    !
    !! if(iOx/=0) call Basis_RedoxAssign(vCpn,vEle,vSpc)
  end select
  !
  call Basis_Build(LBuffer,vCpn)
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
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Change"
  !
end subroutine Basis_Change

subroutine Basis_Change_Wrk(isW,vCpn)
  !---------------------------------------------------------------------
  use M_T_Component,only: T_Component
  use M_Basis_Tools,only: Basis_Cpn_ISpc_Select
  !
  use M_Global_Vars,only: vSpc
  !---------------------------------------------------------------------
  integer,          intent(in)   :: isW
  type(T_Component),intent(inout):: vCpn(:)
  !
  integer,allocatable:: vIndex(:)
  !
  if(iDebug==4) write(fTrc,'(/,A)') "< Basis_Change_wrk"
  !
  allocate(vIndex(size(vCpn)))
  !
  call Basis_Cpn_ISpc_Select( & !
  ! call Basis_Species_Select( & !
  & vSpc,isW,vCpn, & !in
  & vIndex)
  !
  !--------- change basis only if the primary species set has changed --
  if(ANY (vCpn(:)%iSpc /= vIndex(:)) ) then
    !
    vCpn(:)%iSpc= vIndex(:)
    !
    call Basis_Build(.false.,vCpn)
    !
    if(iDebug>2) write(*,'(A)') "===============< basis changed >===="
    !
  end if
  !
  deallocate(vIndex)
  !
  if(iDebug==4) write(fTrc,'(A,/)') "</ Basis_Change_wrk"
end subroutine Basis_Change_Wrk

subroutine Basis_IndexInit(isW,isH_,isOH,isO2,MWsv)
!--
!-- calc' dimensions and indexes that are constant for given vSpc
!--
  use M_T_Species,  only: Species_Index
  !
  use M_Global_Vars,only: vSpc,nAq,nMn,nGs
  !
  integer, intent(out):: isW,isH_,isOH,isO2
  real(dp),intent(out):: MWsv
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_IndexInit"
  !
  nAq= count(vSpc%Typ(1:3)=="AQU")
  nMn= count(vSpc%Typ(1:3)=="MIN")
  nGs= count(vSpc%Typ(1:3)=="GAS")
  nMn= nMn + nGs
  !
  ! find indexes of solvent, H+, OH-, O2 in vSpc
  ! -> unchanged in basis switch -> called only once for a given vSpc
  !-------------- ranks of species H+ and OH- in current species list --
  isW=  Species_Index("H2O",  vSpc)
  isOH= Species_Index("OH-",  vSpc); if(isOH==0) isOH=Species_Index("OH[-]",vSpc)
  isH_= Species_Index("H+",   vSpc); if(isH_==0) isH_=Species_Index("H[+]",vSpc)
  isO2= Species_Index("O2,AQ",vSpc); if(isO2==0) isO2=Species_Index("O2(AQ)",vSpc)
  isO2= Species_Index("O2_AQ",vSpc)
  !
  if(isW==0)  call Stop_("species H2O not found !!!") !--===========stop
  if(isH_==0) call Stop_("species H+ not found  !!!") !--===========stop
  if(isOH==0) call Stop_("species OH- not found !!!") !--===========stop
  !
  MWsv=vSpc(isW)%WeitKg
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_IndexInit"
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
  !
  use M_Global_Vars,only: vSpc
  !
  integer, intent(in)   :: iBal
  integer, intent(in)   :: vOrdPr(:),vOrdAs(:),vOrdMs(:)
  real(dp),intent(inout):: tAlfSp(:,:),tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:),tAlfFs(:,:)
  !
  tAlfSp(iBal,:)=                    real(vSpc(:)%Z)
  tAlfPr(iBal,:)=                    real(vSpc(vOrdPr(:))%Z)
  if(size(vOrdAs)>0) tAlfAs(iBal,:)= real(vSpc(vOrdAs(:))%Z)
  if(size(vOrdMs)>0) tAlfMs(iBal,:)= real(vSpc(vOrdMs(:))%Z)
  tAlfFs(iBal,:)= Zero
  !
end subroutine Basis_InitBalance

subroutine Basis_CpnInert(isW,tStoikio,vCpnIn,vCpnInert)
!--
!-- used only in Equil_Write,
!-- for printing composition of aqueous phase as closed system
!-- make "INERT" all "MOBILE" components,
!-- reorder components, assign prim'species, ...
!--
  use M_T_Component,only: T_Component
  use M_Basis_Tools,only: Basis_Species_Select_Test
  use M_Basis_Tools,only: Basis_Cpn_ISpc_Select
  !
  use M_Global_Vars,only: vSpc
  !
  integer,          intent(in) :: isW
  real(dp),         intent(in) :: tStoikio(:,:)
  type(T_Component),intent(in) :: vCpnIn(:)
  type(T_Component),intent(out):: vCpnInert(:)
  !
  integer:: i
  integer,allocatable:: vIndex(:)
  !
  vCpnInert= vCpnIn
  vCpnInert(:)%Statut="INERT"
  !
  !-------- compute total amount of each component in the fluid phase --
  do i=1,size(vCpnInert)
    vCpnInert(i)%Mole= &
    & SUM( tStoikio(i,:)*vSpc(:)%Dat%Mole, MASK= (vSpc(:)%Typ=="AQU") )
  end do
  !
  !---------------------------------------------- assign prim'species --
  allocate(vIndex(size(vCpnInert)))
  !
  ! call Basis_Species_Select_Test( & !
  ! & vSpc,isW,vCpnInert,tFormula, & !in
  ! & vIndex)

  call Basis_Cpn_ISpc_Select( & !
  ! call Basis_Species_Select( & !
  & vSpc,isW,vCpnInert, & !in
  & vIndex) !out
  !
  vCpnInert(:)%iSpc= vIndex(:)
  !---------------------------------------------/ assign prim'species --
  !
  deallocate(vIndex)
  !
end subroutine Basis_CpnInert

end module M_Basis


