module M_System_Tools
  use M_Kinds
  use M_Trace,only: fTrc,T_,Stop_,iDebug,Pause_,fHtm
  !
  implicit none
  !
  private
  !
  public:: System_Build
  public:: System_Build_Custom
  public:: System_TP_Update
  !
contains

subroutine System_Build
!--
!-- build "master" system --
!-- system with all components (elements) involved in any system in the run
!--
  !!--tools--!!
  use M_T_Species,    only: T_Species,Species_Stoikio_Calc,Species_Append
  use M_Global_Alloc, only: MixModels_Alloc,MixPhases_Alloc,Phases_Alloc !_New
  use M_SolModel_Alloc
  use M_T_SolModel,   only: SolModel_Spc_Init
  use M_T_SolPhase,   only: SolPhase_Init
  use M_Solmodel_Read, only: Solmodel_Solvent_Read
  use M_SolModel_Tools,only: SolModel_TP_Update
  use M_Solmodel_Pitzer_Dtb,only: Solmodel_Pitzer_Dtb_Init,Solmodel_Pitzer_Dtb_TPtest
  use M_Global_Alloc, only: DiscretParam_Alloc
  use M_DiscretModel_Read
  use M_DiscretModel_Tools
  !
  !--global variables--
  use M_Global_Vars,  only: vEle,vSpc,vMixFas,tFormula
  use M_Global_Vars,  only: vMixModel,vDiscretModel,vDiscretParam
  use M_Global_Vars,  only: SolModel
  use M_Global_Vars,  only: vSolModel,vSolFas,vFas
  use M_Global_Vars,  only: nAq,nMn,nGs
  !
  !--database variables--
  use M_Dtb_Vars,      only: DtbFormat,DtbLogK_Dim,DtbLogK_vTPCond
  use M_Solmodel_Vars, only: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  use M_Solmodel_Vars, only: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl,T_Spline
  !
  !--system variables--
  use M_System_Vars,  only: vCpn,TdgK,Pbar
  use M_System_Vars,  only: CpnIsSpc
  use M_System_Vars,  only: System_Zero,System_Type
  !
  integer:: I,N
  logical:: Ok,fOk
  character(len=80):: Msg
  real(dp),allocatable:: vTdgC(:)
  type(T_Species),allocatable:: vSpcTmp(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< System_Build"
  !
  call System_Zero
  !
  !-- read local element list (-current system) -> build vCpn
  call Components_Alloc( & !
  & "SYSTEM",            & !
  & vEle,vSpc,vMixFas,   & !in
  & TdgK,Pbar,           & !
  & System_Type)           !out
  !
  !-- rebuild vEle, with only elements pointed to by vCpn(1:N)
  call Elements_Alloc_forSystem(vCpn)
  !
  ! print *,"debuggg System_TP_Check"
  ! print *,TdgK
  call System_TP_Check(TdgK,Pbar,Ok,Msg)
  if(.not. Ok) call Stop_(trim(Msg))
  !
  if(iDebug>0) then
    do i=1,size(vCpn)
      write(fTrc,'(2A,2(G15.6,1X))') &
      & vCpn(i)%NamCp,vCpn(i)%Statut,vCpn(i)%Mole,vCpn(i)%LnAct
    end do
  end if
  !
  !-- rebuild vSpc, updating vCpn%iSpc
  call Species_Alloc_forSystem(vEle,vCpn)
  !
  !--------------------------------------------- update vSpc(:)%vStoikio
  call Species_Stoikio_Calc(vEle,vSpc,fOk)
  !--------------------------------------------/ update vSpc(:)%vStoikio
  !
  !-------------------------- read solution models -> allocate vMixModel
  call MixModels_Alloc(vSpc,    vMixModel) !
  !
  !------------------------------------------------ add discrete species
  !--- allocate & read vDiscretModel--
  call DiscretModel_Read(vMixModel)
  !
  if(size(vDiscretModel)>0) then
    !
    call DiscretParam_Alloc(vDiscretModel,     vDiscretParam)
    !
    allocate(vSpcTmp(size(vDiscretParam)))
    !
    call DiscretParam_Init( &
    & vEle,vSpc,vMixModel,vDiscretModel, &
    & vDiscretParam,vSpcTmp) !-> build vSpcDiscret
    !
    call Species_Append(vSpcTmp,vSpc) !-> new vSpc !!!
    !
    deallocate(vSpcTmp)
    !
    call DiscretSpecies_Stoikio_Calc( & !
    & vEle,          & !IN
    & vMixModel,     & !IN
    & vDiscretModel, & !IN
    & vDiscretParam, & !IN
    & vSpc)            !INOUT
    !
  end if
  !-----------------------------------------------/ add discrete species
  !
  nAq= count(vSpc%Typ=="AQU")
  nMn= count(vSpc%Typ=="MIN")
  nGs= count(vSpc%Typ=="GAS")
  nMn= nMn + nGs
  !
  !! !--- NEW2010-10-29
  !! N= size(vCpn)
  !! if(CpnIsSpc) then
  !!   do i=1,N
  !!     vCpn(i)%iEle= I
  !!     J= vCpn(i)%iSpc
  !!     vCpn(i)%vStoikCp(0:N+1)= vSpc(J)%vStoikio(0:N+1)
  !!   end do
  !! else
  !!   do i=1,N
  !!     vCpn(i)%iEle= I
  !!     vCpn(i)%vStoikCp(0)= 1 !formula divider !!
  !!     vCpn(i)%vStoikCp(:)= 0
  !!     vCpn(i)%vStoikCp(i)= 1
  !!   end do
  !! end if
  !! !---/ NEW2010-10-29
  !
  !------------------------------------------------ update formula table
  deallocate(tFormula)
  allocate(tFormula(1:size(vEle),1:size(vSpc)))
  call FormulaTable_Calc(vSpc,tFormula)
  !
  if(iDebug>0) call FormulaTable_Sho(vEle,vSpc,tFormula)
  !-----------------------------------------------/ update formula table
  !
  !--- read phase compositions, allocate vMixFas --
  call MixPhases_Alloc(vSpc,vMixModel,  vMixFas)
  call MixPhase_CheckFound(vEle,vSpc,vMixModel,vMixFas,vCpn)
  !
  !----------------------------------------initialize SolModel, SolPhase
  if(System_Type=="AQUEOUS") then
    !
    N= 0
    if(DtbFormat=="LOGKTBL") N= DtbLogK_Dim
    !
    allocate(vTdgC(N))
    if(N>0) vTdgC(1:N)= DtbLogK_vTPCond(1:N)%TdgC
    !
    ! even for HSV base, Solmodel_Read must be called,
    ! for reading activity model
    call Solmodel_Solvent_Read( &
    & N,vTdgC,SolModel, &
    & Ok_Rho, Ok_Eps, Ok_DHA, Ok_DHB, Ok_BDot, &
    & Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl)
    !
    deallocate(vTdgC)
    !
    if(SolModel%iActModel== 8  .or. &    !! "PITZER"
    &  SolModel%iActModel== 11 ) then    !! "SIT"
      call Solmodel_Pitzer_Dtb_Init(vSpc)
    end if
    !
    call SolModel_TP_Update(TdgK,Pbar,SolModel)

    if(iDebug>2) then
      if(SolModel%iActModel== 8  .or. &    !! "PITZER"
      &  SolModel%iActModel== 11 ) &    !! "SIT"
      & call Solmodel_Pitzer_Dtb_TPtest !!!"WRK"!!!
    end if
    !
    !----------------------------------------------------------------NEW
    !--- initialize indexes of the sol'model
    call SolModel_Spc_Init("H2O",vSpc,SolModel,Ok,Msg)
    call SolModel_Alloc ! allocate vSolModel
    vSolModel(1)= SolModel
    !
    call SolPhase_Alloc ! allocate vSolFas
    call SolPhase_Init( & !
    & 1,          & !IN, model index
    & vSolModel,  & !IN, base of solution models
    & vSolFas(1), & !INOUT, solution phase
    & Ok,Msg)       !OUT
    !
    !! call SolPhase_Model_Init(SolModel%Model,vSolModel,vSolFas(1),Ok,Msg)
    !! if(.not. Ok) call Stop_(trim(Msg))
    !
    if(.not. Ok) &
    & call Stop_("SolModel_Spc_Init "//trim(Msg))
    !---------------------------------------------------------------/NEW

  end if
  !---------------------------------------/initialize SolModel, SolPhase
  !
  !-- allocate vFas
  ! call Phases_Alloc_New(vSpc,vMixFas,vSolFas)
  call Phases_Alloc(vSpc,vMixFas,  vFas)
  !
  !------------------------------------- update global thermo'parameters
  !--------------------- (vSpc,vMixModel,vMixFas,vFas) at (TdgK,Pbar) --
  call System_TP_Update(TdgK,Pbar)
  !------------------------------------/ update global thermo'parameters
  !
  if(System_Type=="AQUEOUS") then
    !-- check presence of species H2O, H+, OH- in vSpc --
    call System_Species_Check(vSpc,Ok,Msg)
    if(.not. Ok) call Stop_(Msg)
  end if
  !
  !! call Basis_Init(vCpn,TdgK,Pbar)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ System_Build"
  !
end subroutine System_Build

subroutine Species_Alloc_forSystem(vEle,vCpn)
!--
!-- from general vSpc built from database,
!-- reduce vSpc to species involved in the run
!--
  use M_Numeric_Const,       only: Ln10
  use M_T_Component, only: T_Component
  use M_T_Element,   only: T_Element,Element_Index,Formula_Read
  use M_T_Species,   only: T_Species,Species_Index
  !
  use M_Global_Vars,only: vSpc,tFormula
  use M_System_Vars,only: System_Type
  !
  type(T_Element),   intent(inout):: vEle(:)
  type(T_Component), intent(inout):: vCpn(:)
  !
  character(len=23),dimension(:),allocatable:: vNamSpc
  !
  integer,allocatable:: vStoik(:)
  !! logical,allocatable:: vExclude(:)
  !! logical,allocatable:: vInclude(:)
  !
  type(T_Species),allocatable:: vSpcNew(:)
  !
  type(T_Species):: S_
  logical        :: fOk
  integer        :: nCp,iCp,nSp,iSp,iOx_,ZSp,Z_,i,nAq,nMn,nGs
  !
  if(iDebug>0) write(fTrc,'(/,A,/)') "< Species_Alloc"
  !
  nCp=size(vCpn)
  !
  !before restructuring vSpc, save vCpn-vSpc links
  allocate(vNamSpc(1:nCp))
  do iCp=1,nCp
    vNamSpc(iCp)= trim(vSpc(vCpn(iCp)%iSpc)%NamSp)
  end do
  !
  call Redox_Calc(vSpc, vEle,vCpn)
  !
  !! !------------------ read SPECIES.include and SPECIES.EXCLUDE blocks --
  !! allocate(vExclude(1:size(vSpc)))  ;  vExclude(:)= .false.
  !! allocate(vInclude(1:size(vSpc)))  ;  vInclude(:)= .true.
  !! call Species_Read_Excluded(vSpc,vExclude,vInclude)
  !! !-----------------/ read SPECIES.include and SPECIES.EXCLUDE blocks --
  !
  !----------------------------------------------------species selection
  !-- from the vSpc built from database,
  !-- retrieve species consistent
  !-- with vEle & redox
  if(iDebug>0) write(fTrc,'(A)') &
  & "< Restrict database to species consistent with element list and redox state"
  !
  allocate(vSpcNew(size(vSpc)))
  !
  allocate(vStoik(1:nCp))
  iOx_= Element_Index("OX_",vEle)
  nSp=0
  do iSp=1,size(vSpc)
    !
    !! if( vExclude(iSp) .or. (.not. vInclude(iSp)) ) cycle
    !
    S_=vSpc(iSp)
    if(S_%NamSp(1:1)=="!") cycle !!!§§§ bricollage !!!
    !
    call Formula_Read(S_%Formula,vEle,ZSp,S_%Div,fOk,vStoik)
    !
    if(fOk) then
    !fOk=FormulaIsOk (as regards elements in run), build new LnkSpc, called Lnk
      !
      Z_=dot_product(vStoik(1:nCp),vEle(1:nCp)%Z) !check charge of species
      !
      if(iOx_/=0 .or. ZSp==Z_) then !if charge is OK, then save in vSpcNew
        !
        nSp=  nSp+1
        S_%Z= ZSp
        !
        vSpcNew(nSp)= S_
        !
        if(iDebug>0) write(fTrc,'(A,A1,A15,A1,A39,A1,2I3)') &
        & "ACCEPT",T_,S_%NamSp,T_,S_%Formula,T_,S_%Z,nSp
      end if
      !
    end if
    !
  end do
  !
  !! deallocate(vExclude)
  !! deallocate(vInclude)
  !! deallocate(vStoik)
  !--------------------------------------------------/ species selection
  !
  !----------------------------------------------- build new sorted vSpc
  deallocate(vSpc)
  allocate(vSpc(1:nSp))
  !
  nAq= 0; nMn= 0; nGs= 0
  iSp= 0
  !--retrieve first aqu'species
  do i=1,nSp
    if(vSpcNew(i)%Typ=="AQU") then
      iSp= iSp +1
      vSpc(iSp)= vSpcNew(i)
    end if
  end do
  !--then retrieve first min'species
  do i=1,nSp
    if(vSpcNew(i)%Typ=="MIN") then
      iSp= iSp +1
      vSpc(iSp)= vSpcNew(i)
    end if
  end do
  !--then retrieve first gas'species
  do i=1,nSp
    if(vSpcNew(i)%Typ=="GAS") then
      iSp= iSp +1
      vSpc(iSp)= vSpcNew(i)
    end if
  end do
  deallocate(vSpcNew)
  !----------------------------------------------/ build new sorted vSpc
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Restrict database"
  !
  if(count(vSpc%Typ=="AQU")==0 .and. System_Type=="AQUEOUS") &
  & call Stop_("Found NO Aqueous Species ...") !--------------------stop
  !
  !-------------------------------------------then relate vCpn with vSpc
  if(iDebug>0) write(fTrc,'(/,A,/)') &
  & "Species_Alloc_forSystem, Element/ Component/ Status/ iSpc"
  !
  !vSpc may have changed
  !-> need to check whether the species pointed to by Cpn's are in the new vSpc
  !
  do iCp=1,nCp
  !---------------------find which vSpc(iSp) matches with vNamSpc(iCp)--
    !
    iSp= Species_Index(vNamSpc(iCp),vSpc) !-> index of Prim'species in new vSpc
    !
    if(iSp==0) call Stop_(trim(vNamSpc(iCp))//" Not Found as Species")
    !
    if(System_Type=="AQUEOUS") then
      if(vSpc(iSp)%Typ /= "AQU") then
        if(trim(vCpn(iCp)%Statut)=="INERT") &
        & call Stop_( & !-------------------------------------------stop
        & trim(vNamSpc(iCp))//" : ALL INERT SPECIES SHOULD BE AQUEOUS")
      end if
    end if
    !
    !update %iSpc (link between vSpc and vCpn)
    vCpn(iCp)%iSpc= iSp
    !
    if(iDebug>0) write(fTrc,'(A3,A23,A7)') &
    & vEle(iCp)%NamEl,vNamSpc(iCp),vCpn(iCp)%Statut
    !
  end do
  !
  deallocate(vNamSpc)
  !
  if(iDebug>1) then
    write(fTrc,'(/,A,/)') "Species list, final"
    do iSp=1,nSp
      write(fTrc,'(I7,A1,A23,A1,A3)') iSp,T_,vSpc(iSp)%NamSp,T_,vSpc(iSp)%Typ
    end do
  end if
  !
  if(iDebug>0) write(fTrc,'(/,A,/)') "</ Species_Alloc"
  !
end subroutine Species_Alloc_forSystem

subroutine System_Species_Check(vSpc,Ok,Msg)
!--
!-- check presence of species H2O, H+, OH- in vSpc --
!-- (for an aqueous system)
!--
  use M_T_Species,only: T_Species,Species_Index
  !
  type(T_Species),  intent(in)  :: vSpc(:)
  logical,          intent(out) :: Ok
  character(len=80),intent(out):: Msg
  !
  Ok=  .true.
  Msg= "Ok"
  !
  if(Species_Index("H2O",vSpc)==0) then
    Msg= "species H2O not found !!!"
    Ok=  .false.
    return
  end if
  if(Species_Index("OH-",vSpc)==0 .and. Species_Index("OH[-]",vSpc)==0) then
    Msg= "species OH- not found !!!"
    Ok=  .false.
    return
  end if
  if(Species_Index("H+", vSpc)==0 .and. Species_Index("H[+]",vSpc)==0) then
    Msg= "species H+ not found !!!"
    Ok=  .false.
    return
  end if
  !
end subroutine System_Species_Check

subroutine System_TP_Check(TdgK,Pbar,Ok,Msg)
  use M_T_Tpcond,only: T_TPCond
  use M_Dtb_Vars,only: DtbFormat,DtbLogK_vTPCond,Psat_Auto
  use M_Dtb_Calc,only: Dtb_TP_Check
  !
  real(dp),        intent(in)   :: TdgK
  real(dp),        intent(inout):: Pbar
  logical,         intent(out)  :: Ok
  character(len=*),intent(out)  :: Msg
  !
  Ok=  .true.
  Msg= "OK"
  call Dtb_TP_Check(DtbFormat,DtbLogK_vTPCond,Psat_Auto,TdgK,Pbar,Ok,Msg)
  !
end subroutine System_TP_Check

subroutine System_TP_Update(TdgK,Pbar)
!--
!-- update thermo parameters --
!--
  use M_Global_Tools,  only: Global_TP_Update
  use M_SolModel_Tools,only: SolModel_TP_Update
  !
  use M_Global_Vars, only: vSpcDtb,vSpc,vMixModel,vDiscretModel,vDiscretParam
  use M_Global_Vars, only: vMixFas,vFas,SolModel
  use M_System_Vars, only: System_Type
  !
  real(dp),intent(in):: TdgK,Pbar

  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)

  if(System_Type=="AQUEOUS") call SolModel_TP_Update(TdgK,Pbar,SolModel)

end subroutine System_TP_Update

subroutine MixPhase_CheckFound( &
& vEle,vSpc,vMixModel, &
& vMixFas,vCpn)
!--
!-- check whether phases have been found
!-- for all components with %namSol /- "AQU"
!--
  use M_T_Component,only: T_Component
  use M_T_Element,  only: T_Element
  use M_T_Species,  only: T_Species
  use M_T_MixModel, only: T_MixModel,Mix_Site
  use M_T_MixModel, only: MixModel_Site_ActivIdeal,MixModel_XPoleToXSite
  use M_T_MixPhase, only: T_MixPhase, MixPhase_Index
  !
  type(T_Element),  intent(in)   :: vEle(:)
  type(T_Species),  intent(in)   :: vSpc(:)
  type(T_MixModel), intent(in)   :: vMixModel(:)
  type(T_MixPhase), intent(inout):: vMixFas(:) !-> update vMixFas%vLPole
  type(T_Component),intent(inout):: vCpn(:) !-> update vCpn%iPol, vCpn%iMix
  !
  type(T_MixModel):: SM
  integer :: nCp,iSp,iCp,iMix,J,K
  real(dp):: Act
  LOGICAl, allocatable:: Phase_Found(:)
  real(dp),allocatable:: vXAtom(:)
  !-> check that phase for each component is described in input
  !
  if(iDebug>0) write(fTrc,'(/,A,/)') "< MixPhases_CheckFound"
  !
  nCp= size(vCpn)
  !
  allocate(Phase_Found(1:nCp))  ;  Phase_Found=.false.
  !
  do iCp=1,nCp !check components
    !
    iSp=vCpn(iCp)%iSpc
    !
    if(iDebug>0) write(fTrc,'(I3,6A)') &
    & iCp, &
    & " >> ELE=",  vEle(vCpn(iCp)%iEle)%NamEl, &
    & " >> SPC=",  vSpc(vCpn(iCp)%iSpc)%NamSp, &
    & " >> PHASE=",trim(vCpn(iCp)%namSol)
    !
    if(iSp==0) call Stop_ &
    & ("NO SPECIES FOR COMPONENT "//trim(vEle(vCpn(iCp)%iEle)%NamEl))
    !
    if(trim(vCpn(iCp)%namSol)=="Z") then
      !
      Phase_Found(iCp)=.true.
      !
    else
      !
      iMix= MixPhase_Index(vMixFas,vCpn(iCp)%namSol)
      !
      if(iMix==0) call Stop_ & !------------------------------------stop
      & ("PHASE NOT FOUND FOR COMPONENT "//trim(vEle(vCpn(iCp)%iEle)%NamEl))
      !
      Phase_Found(iCp)=.true.
      vCpn(iCp)%iMix= iMix !!vSpc(iSp)%iMix=iMix
      !
      if(iDebug>0) write(fTrc,'(4(A,1X))') &
      & "COMPONENT", vEle(vCpn(iCp)%iEle)%NamEl, &
      & "PHASE", trim(vMixFas(vCpn(iCp)%iMix)%Name)
      !
      !---------------------find index of "species" in "solution"%vIPole
      SM= vMixModel(vMixFas(iMix)%iModel)
      J=0
      do K=1,SM%NPole
        if(SM%vIPole(K)==iSp) J=K
      end do
      !-----------------/
      !
      if(iDebug>0) write(fTrc,'(2A,/)') &
      & "END-MEMB  ", trim(vSpc(iSp)%NamSp)
      !
      if(J==0) call Stop_& !----------------------------------------stop
      & ("SPECIES NOT IN MIXTURE:"//trim(vSpc(iSp)%NamSp))
      !
      vCpn(iCp)%iPol=J
      !--- iPol- index of species in end-member list of phase
      if(SM%Model==Mix_Site) then
        allocate(vXAtom(SM%NAtom))
        call MixModel_XPoleToXSite( &
        & SM, vMixFas(iMix)%vXPole, &
        & vXAtom)
        !
        Act= MixModel_Site_ActivIdeal(SM,J,vXAtom)
        !
        if(iDebug>0) &
        & write(fTrc,'(A,A,G15.6)') trim(vSpc(iSp)%NamSp),">> ACT=",Act
        !
        if(Act/=0.0D0) vMixFas(iMix)%vLPole=.true.
        !
        deallocate(vXAtom)
      end if
      !
    end if
    !
  end do
  !
  do J=1,nCp
    if(.not. Phase_Found(J)) & !------------------------------------stop
    & call Stop_( "DESCRIPTION NOT FOUND FOR "//trim(vCpn(J)%namSol) )
  end do
  !
  deallocate(Phase_Found)
  !
  if(iDebug>0) write(fTrc,'(/,A,/)') "</ MixPhases_CheckFound"
end subroutine MixPhase_CheckFound

subroutine Components_Alloc(  & !
!--
!-- read SYSTEM block, build vCpn
!-- (using temporary vCpn)
!--
& sKeyWord,vEle,vSpc,vMixFas, & !IN
& TdgK,Pbar,                  & !INOUT
& SysType)                      !OUT
  use M_T_Element,  only: T_Element
  use M_T_Species,  only: T_Species
  use M_T_Component,only: T_Component,Component_Zero
  use M_T_Tpcond,   only: T_TPcond
  use M_T_MixPhase, only: T_MixPhase
  use M_Component_Read,only: Components_Read
  !
  use M_System_Vars,only: vCpn
  !
  character(len=*),intent(in):: sKeyWord
  type(T_Element), intent(in):: vEle(:)
  type(T_Species), intent(in):: vSpc(:)
  type(T_MixPhase),intent(in):: vMixFas(:)
  !
  real(dp),        intent(inout):: TdgK,Pbar
  character(len=7),intent(out)  :: SysType
  !
  integer,parameter:: MaxCpn= 60  ! dimension of temporary vCpn
  !
  type(T_Component):: vTmp(MaxCpn)
  character(len=80):: Msg
  logical:: Ok
  integer:: N,I
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Components_Alloc"
  !
  do I=1,size(vTmp)
    call Component_Zero(vTmp(I))
  end do

  call Components_Read( & !
  & sKeyWord,   & !IN
  & vEle,vSpc,vMixFas, & !IN
  & TdgK,Pbar,  & !INOUT
  & N,          & !OUT
  & SysType,    & !OUT
  & Ok,         & !OUT
  & Msg,        & !OUT
  & vTmp)         !OUT
  !
  if(.not.Ok) call Stop_(trim(Msg))
  if(N==0) call Stop_("NO COMPONENTS FOUND")

  if(allocated(vCpn)) deallocate(vCpn)
  allocate(vCpn(1:N))
  vCpn(1:N)= vTmp(1:N)

  if(iDebug>0) then
    do i=1,N
      write(fTrc,'(A)') vCpn(I)%NamCp
    end do
  end if

  if(iDebug>0) write(fTrc,'(A,/)') "</ Components_Alloc"

  return
end subroutine Components_Alloc

subroutine System_Build_Custom( &
!--
!-- from block sKeyWord, build a (sub)system vCpn of master system vCpnMaster
!--
& sKeyWord,          & !IN
& vEle,vSpc,vMixFas, & !IN
& vCpnMaster,        & !IN
& TdgK,Pbar,         & !INOUT
& vCpn,              & !INOUT
& Ok)                  !OUT
!--
!-- CAVEAT=
!-- the input master system must have all components INERT (or BUFFER) !!WHY ??
!--
  use M_T_Element,  only: T_Element
  use M_T_Species,  only: T_Species
  use M_T_Tpcond,   only: T_TPcond
  use M_T_MixPhase, only: T_MixPhase
  use M_T_Component,only: T_Component,Component_Index,Component_Print,CpnMolMinim
  use M_Component_Read,only: Components_Read
  !
  character(len=*), intent(in):: sKeyWord
  type(T_Element),  intent(in):: vEle(:)
  type(T_Species),  intent(in):: vSpc(:)
  type(T_MixPhase), intent(in):: vMixFas(:)
  type(T_Component),intent(in):: vCpnMaster(:)
  !
  real(dp),         intent(inout):: TdgK,Pbar
  type(T_Component),intent(inout):: vCpn(:)
  logical,          intent(out)  :: Ok
  !
  character(len=7):: SysType
  !
  logical,          allocatable:: vCpnFound(:)
  type(T_Component),allocatable:: vCpnNew(:)
  !
  type(T_Component):: C
  integer:: iCp,I,N,nCp
  !!logical:: Ok
  character(len=80):: Msg
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< System_Build_Custom"
  !
  !-------------------default values for component not found in BuildLnk
  vCpn= vCpnMaster  !! <-must be done "outside" ???
  vCpn(:)%Mole=   CpnMolMinim
  where(vCpn(:)%Statut/= "BUFFER") vCpn(:)%Statut= "INERT"
  !--------------------------------------------------------------------/
  !
  nCp= size(vCpn)
  allocate(vCpnFound(1:nCp))  ;  vCpnFound=  .false.
  allocate(vCpnNew(1:nCp))    ;  vCpnNew(:)= vCpn(:)
  !
  call Components_Read( & !
  & sKeyWord,   & !IN
  & vEle,vSpc,vMixFas, & !IN
  & TdgK,Pbar,  & !INOUT
  & N,          & !OUT
  & SysType,    & !OUT
  & Ok,         & !OUT
  & Msg,        & !OUT
  & vCpnNew)      !OUT
  !
  !~ if(.not. Ok) call Stop_(trim(Msg))
  if(.not. Ok) return
  !
  
  ! print *,"debuggg System_Build_Custom"
  ! print *,TdgK
  call System_TP_Check(TdgK,Pbar,Ok,Msg)
  !~ if(.not. Ok) call Stop_(trim(Msg))
  if(.not. Ok) return
  !
  !if(N==0) call Stop_("NO COMPONENTS FOUND")
  Ok= N>0
  if(.not. Ok) return !------------------------------------N=0 -> return

  !--------------------------------------components listed in new system
  do I=1,N
    C= vCpnNew(I)
    iCp= Component_Index(C%NamCp,vCpnMaster)
    if(iCp==0) then
      Ok= .false.
      call Stop_("Error in System_Build_Custom: new Cpn not in vCpnMaster")
    end if
    !!print *, "Ele,iCp=",vEle(C%iEle)%Name, iCp  ;  pause
    !-> will have same element order in vCpn as in vCpnMaster
    vCpnFound(iCp)=.true.
    !! if(C%Statut=="BUFFER") then
    !!   if(vCpnMaster(iCp)%Statut/="BUFFER") &
    !!   & call Stop_ &
    !!   & ("A Cpn cannot be BUFFER in subsystem if it's not also BUFFER in master system !!!")
    !! end if
    !if(vCpnMaster(iCp)%Statut/="BUFFER")
    vCpn(iCp)= C
    !-> when a Cpn is BUFFER in master system, it is BUFFER in any system
  end do

  !--------------------------------- components not listed in new system
  !old- when an element is not listed, it is "absent": %Mole-CpnMolMinim
  !new- when an element is not listed, it is like in the master system
  do iCp=1,nCp
    if(.not. vCpnFound(iCp)) then
      !
      vCpn(iCp)= vCpnMaster(iCp)
      !
      if(vCpnMaster(iCp)%Statut/="BUFFER") then
        !
        vCpn(iCp)%Statut= "INERT"
        !
        if(vSpc(vCpn(iCp)%iSpc)%Typ /= "AQU") then
          print '(A)', "unlisted components in subsystem are INERT !!"
          print '(2A)',vSpc(vCpn(iCp)%iSpc)%NamSp,"-> species should be aqueous !!!"
          Ok= .false.
        end if
      end if
      !
    end if
  end do
  !
  deallocate(vCpnNew)
  deallocate(vCpnFound)
  !
  !--------update global thermo'parameters (vSpc,vMixModel,vMixFas,vFas)
  !--------at (TdgK,Pbar)
  call System_TP_Update(TdgK,Pbar)
  !-------------------------------------/update global thermo'parameters
  
  !----------------------------------------------------------------trace
  if(iDebug>2) then
    print '(A)',"< System_Build_Custom "//trim(sKeyWord)
    do I=1,size(vCpn)
      call Component_Print(6,vEle,vSpc,vCpn(I))
    end do
    print '(A)',"</ System_Build_Custom"
    call Pause_
  end if
  !---------------------------------------------------------------/trace
  
  if(iDebug>0) write(fTrc,'(A,/)') "</ System_Build_Custom"
  
end subroutine System_Build_Custom

subroutine Elements_Alloc_forSystem(vCpn)
!--
!-- rebuild vEle, with only elements pointed to by vCpn(1:N)
!--
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  use M_Global_Vars,only: vEle
  !
  type(T_Component),intent(inout):: vCpn(:)
  !
  type(T_Element),allocatable:: vEleB(:)
  integer:: I,N
  !
  N= size(vCpn)
  !
  allocate(vEleB(1:size(vEle)))
  vEleB= vEle !-------------------------------save current vEle to vEleB
  !
  deallocate(vEle) !----------------------------------- destroy old vEle
  allocate(vEle(1:N))
  vEle(1:N)= vEleB(vCpn(1:N)%iEle) !----------------------build new vEle
  !
  !----------------------------- update fields in component's constructs
  do i=1,N
    vCpn(i)%iEle= I
    vCpn(i)%vStoikCp(0)= 1 !formula divider !!
    vCpn(i)%vStoikCp(:)= 0
    vCpn(i)%vStoikCp(i)= 1
  end do
  !--------------------------------------------------------------------/
  !
  deallocate(vEleB)
  !
end subroutine Elements_Alloc_forSystem

subroutine Redox_Calc(vSpc,vEle,vCpn)
!--
!-- first, from the names of the primary species,
!-- compute the valencies of "potentially redox elements" (Fe,Cl,...)
!--
  use M_T_Component, only: T_Component
  use M_T_Element,   only: T_Element,Element_Index,Formula_Read
  use M_T_Species,   only: T_Species,Species_Stoikio
  !
  type(T_Species),   intent(in)   :: vSpc(:)
  type(T_Element),   intent(inout):: vEle(:)
  type(T_Component), intent(inout):: vCpn(:)
  !
  type(T_Species):: S
  integer:: ieOx, iCp, nCp, Z !, Zsp
  logical:: fOk
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Redox_Calc"
  !
  ieOx= Element_Index("OX_",vEle)
  !
  nCp= size(vEle)
  !
  do iCp=1,nCp
    !
    ! restrict this procedure to elements with no valency
    ! assigned in dtb (Fe,S,Cr,As)
    if(vEle(iCp)%Z==0) then !or Ele%ReDox=="VAR" ???

      S= vSpc(vCpn(iCp)%iSpc)

      call Species_Stoikio(vEle,ieOx,S,fOk)
      ! nota: true stoikiometry is vStoik/S%Div, true charge is Zsp/S%Div
      
      ! if(abs(vStoik(iCp))/=1) &
      ! & call Stop_("coeff' should be +/-1")

      ! if(fOk) then !fOk=FormulaIsOk (as regards elements in run)
      ! calculate reDox state of element in its Prim'Species -> basis valency

      Z= dot_product( S%vStoikio(1:nCp),  vEle(1:nCp)%Z ) &
      &             - S%vStoikio(iCp  ) * vEle(iCp  )%Z
      !
      vEle(iCp)%Z= (S%Z-Z) / S%vStoikio(iCp)
      !
      if(iDebug>0) write(fTrc,'(2I3,A,A3,A1,3I3)') &
      & S%Z-Z, S%vStoikio(iCp), &
      & " -> Default Valency of ", vEle(iCp)%NamEl, "=", vEle(iCp)%Z
      
      !end if

    end if
  end do
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Redox_Calc"
  !
end subroutine Redox_Calc
  !
subroutine FormulaTable_Calc(vSpc,T)
!--
!-- compute the formula table elements / species
!--
  use M_T_Species,  only: T_Species
  !
  type(T_Species),intent(in) :: vSpc(:)
  real(dp),       intent(out):: T(:,:)
  !
  integer:: I, N
  !
  if(iDebug>0) write(fTrc,'(/,A,/)') "< FormulaTable_Calc"
  !
  N= size(T,1)
  do I=1,size(vSpc)
    !J= vSpc(I)%vStoikio(0)
    !if(J==0) J=1
    T(1:N,I)= vSpc(i)%vStoikio(1:N) / real(vSpc(I)%vStoikio(0))
    !column I of T = Formula of species I
    !-> T is the "Formula Matrix" a la SmithMissen
  end do
  !
  if(iDebug>0) write(fTrc,'(/,A,/)') "</ FormulaTable_Calc"
  !
end subroutine FormulaTable_Calc

subroutine FormulaTable_Sho(vEle,vSpc,T)
  use M_IoTools,  only: GetUnit
  use M_Files,    only: DirOut,Files_Index_Write
  use M_T_Element,only: T_Element
  use M_T_Species,only: T_Species
  !
  type(T_Element),  dimension(:),  intent(in):: vEle
  type(T_Species),  dimension(:),  intent(in):: vSpc
  real(dp),         dimension(:,:),intent(in):: T
  !
  integer:: F,I,iEl
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_stoik.log")
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_stoik.log",&
  & "STOIKIO: formula matrix tStoikio, transposed")
  !
  do iEl=1,size(vEle)
    write(F,'(A3,A1)',advance="no") vEle(iEl)%NamEl, T_
  end do
  write(F,'(3(A,A1))') "Chg",T_,"species",T_
  !
  do I=1,size(vSpc)
    do iEl=1,size(vEle)
      write(F,'(F7.2,A1)',advance="no") T(iEl,I), T_
    end do
    write(F,'(I3,A1,A12,A1)') vSpc(I)%Z, T_, trim(vSpc(I)%NamSp),T_
  end do
  !
  close(F)
end subroutine FormulaTable_Sho

end module M_System_Tools
