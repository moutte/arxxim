module M_Global_Build
  implicit none
  !
  private
  !
  public:: Global_Build
  ! public:: Global_Build_New

contains

subroutine Global_Build
!--
!-- build "global" vEle, vSpc, vMixModel, vMixFas, vFas
!-- (all available elements and species, not restricted to those used in the run)
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc  !,Stop_,Warning_
  !
  use M_T_Species,   only: T_Species,Species_Stoikio_Calc
  use M_Dtb_Const,   only: T_CK,Tref,Pref
  use M_Dtb_Read,    only: Dtb_Read_Species
  use M_Dtb_Calc,    only: SpeciesDtb_ToSpecies
  use M_Element_Read,only: Element_Read_Redox,Element_Read_Entropy
  use M_SpeciesDtb_Read,only: Species_Read_AquSize,Species_Write_AquSize
  use M_Global_Tools,only: Global_TP_Update,Global_Species_Select
  use M_Global_Alloc
  !
  !--global variables--
  use M_Global_Vars, only: vEle,vSpc,vSpcDtb
  use M_Global_Vars, only: vMixModel,vMixFas,vFas
  use M_Global_Vars, only: vDiscretModel,vDiscretParam
  use M_Global_Vars, only: vSolModel
  use M_Global_Vars, only: MySpace
  !
  integer :: iTP, N
  logical :: fOk
  real(dp):: TdgK,Pbar
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Global_Build"
  !
  call Elements_Alloc_forDtb(vEle,N) !---------------------> builds vEle
  !
  call Element_Read_Redox(vEle)
  !call Element_Read_Entropy(vEle)
  !
  call Dtb_Read_Species(vEle) !--------> read databases -> build vDtbXxx
  !
  call SpeciesDtb_Alloc(vSpcDtb,N) !---------> from vDtb*, build vSpcDtb
  !
  ! write(11,'(A)') "<=======vSpcDtb=1==="
  ! do n=1,size(vSpcDtb)
  !   write(11,'(I3,1X,I3,1X,A)') n, vSpcDtb(n)%Indx, vSpcDtb(n)%DtbModel
  ! end do
  ! write(11,'(A)') "===================>"
  !
  N= size(vSpcDtb)
  if(N<1) return
  !
  deallocate(vSpc)
  allocate(vSpc(N))
  !
  call SpeciesDtb_ToSpecies(vSpcDtb,vSpc)
  !
  ! write(11,'(A)') "<=======species=2==="
  ! do n=1,size(vSpc)
  !   write(11,'(I3,1X,A)') n, vSpc(n)%NamSp
  ! end do
  ! write(11,'(A)') "===================>"
  !
  call Global_Species_Select
  !
  call Species_Stoikio_Calc(vEle,vSpc,fOk)
  !
  call Species_Read_AquSize(vSpc)
  ! call Species_Write_AquSize(vSpc)
  !
  call vSpecies_Rename(vSpc)
  !
  call MixModels_Alloc(vSpc,  vMixModel) !-------build MixModel database
  !
  call MixPhases_Alloc(vSpc,vMixModel,  vMixFas)
  !
  call Phases_Alloc(vSpc,vMixFas,  vFas) !------------------> build vFas
  !
  call Condition_Read(iTP)
  !
  TdgK= Tref
  Pbar= Pref
  !
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
  if(.not.allocated(vSolModel)) allocate(vSolModel(0)) !!JM!!EN_COURS
  !
  if(allocated(MySpace%vEle)) call MySpace%Free_ !!JM!!
  call MySpace%New_( &
  & size(vEle),      & !
  & size(vSpcDtb),   & !
  & size(vSpc),      & !
  & size(vMixModel), & !
  & size(vSolModel)  )
  !!JM!!EN_COURS ... MySpace must be free'd at end of loop ....
  !
  call Species_Write_AquSize(vSpc)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Global_Build"
  !
end subroutine Global_Build

subroutine Global_Build_New !-------------------------------------UNUSED
!--
!-- build "global" vEle, vSpc, vMixModel, vMixFas, vFas
!-- (all available elements and species,
!-- not restricted to those used in the run)
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,Stop_  !,Warning_
  !
  use M_T_Species,   only: T_Species,Species_Stoikio_Calc,Species_Append
  use M_T_SolModel,  only: T_SolModel
  !
  use M_Dtb_Const,   only: T_CK,Tref,Pref
  use M_Dtb_Read,    only: Dtb_Read_Species
  use M_Dtb_Calc,    only: SpeciesDtb_ToSpecies
  !
  use M_Element_Read,only: Element_Read_Redox,Element_Read_Entropy
  use M_SpeciesDtb_Read,only: Species_Read_AquSize,Species_Write_AquSize
  use M_Global_Tools,only: Global_TP_Update,Global_Species_Select
  use M_Global_Alloc
  !
  use M_SolModel_Read
  use M_SolPhase_Read
  use M_T_SolModel,   only: SolModel_Spc_Init
  use M_T_SolPhase,   only: SolPhase_Init
  !
  use M_DiscretModel_Read
  use M_DiscretModel_Tools
  !
  ! !--database variables--
  ! use M_Solmodel_Vars, only: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  ! use M_Solmodel_Vars, only: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl,T_Spline
  !
  !--global variables--
  use M_Global_Vars, only: vEle,vSpc,vSpcDtb
  use M_Global_Vars, only: vMixModel,vMixFas,vFas
  use M_Global_Vars, only: vSolModel,vSolFas
  use M_Global_Vars, only: vDiscretModel,vDiscretParam
  !
  integer :: iTP, N
  logical :: fOk
  real(dp):: TdgK,Pbar
  logical :: Ok
  character(len=80):: Msg
  type(T_Species),allocatable:: vSpcTmp(:)
  ! type(T_SolModel):: SolModel
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Global_Build"
  !
  call Elements_Alloc_forDtb(vEle,N) !---------------------> builds vEle
  !
  call Element_Read_Redox(vEle)
  !call Element_Read_Entropy(vEle)
  !
  call Dtb_Read_Species(vEle) !--------> read databases -> build vDtbXxx
  !
  call SpeciesDtb_Alloc(vSpcDtb,N) !---------> from vDtb*, build vSpcDtb
  !
  N= size(vSpcDtb)
  if(N<1) return
  !
  deallocate(vSpc)
  allocate(vSpc(1:N))
  !
  call SpeciesDtb_ToSpecies(vSpcDtb,vSpc)
  !
  call Global_Species_Select
  !
  call Species_Stoikio_Calc(vEle,vSpc,fOk)
  !
  call Species_Read_AquSize(vSpc)
  !call Species_Write_AquSize(vSpc)
  !
  call vSpecies_Rename(vSpc)
  !
  call MixModels_Alloc(vSpc, vMixModel) !--------build MixModel database
  !
  !-------------------------------------------------add discrete species
  !--allocate & read vDiscretModel--
  call DiscretModel_Read(vMixModel)
  !
  if(size(vDiscretModel)>0) then
    !
    call DiscretParam_Alloc(vDiscretModel,  vDiscretParam)
    !
    allocate(vSpcTmp(size(vDiscretParam)))
    !
    call DiscretParam_Init( &
    & vEle,vSpc,vMixModel,vDiscretModel, &
    & vDiscretParam,vSpcTmp) !-> build vSpcDiscret
    !
    call Species_Append(vSpcTmp,      vSpc) !-> new vSpc !!!
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
  !------------------------------------------------/add discrete species
  !
  !----------------------------------------------------- build vSolModel
  if(count(vSpc(:)%Typ=="AQU") >0) then

    ! N= 0
    ! if(DtbFormat=="LOGKTBL") N= DtbLogK_Dim

    ! allocate(vTdgC(N))
    ! if(N>0) vTdgC(1:N)= DtbLogK_vTPCond(1:N)%TdgC

    ! ! even for HSV base, Solmodel_Read must be called,
    ! ! for reading activity model
    ! call Solmodel_Solvent_Read( &
    ! & N,vTdgC,SolModel, &
    ! & Ok_Rho, Ok_Eps, Ok_DHA, Ok_DHB, Ok_BDot, &
    ! & Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl)

    ! if(SolModel%ActModel=="PITZER" .or. &
    ! &  SolModel%ActModel=="SIT"  ) then
    !   call Pitzer_Dtb_Init(vSpc)
    !   call Pitzer_Dtb_TPUpdate(TdgK,Pbar)
    !   if(iDebug>2) call Pitzer_Dtb_TPtest !!!"WRK"!!!
    ! end if

    !!--- initialize indexes of the sol'model
    ! call SolModel_Spc_Init("H2O",vSpc,SolModel,Ok,Msg)

    ! call SolModel_Alloc ! allocate vSolModel
    ! vSolModel(1)= SolModel

    call SolModel_Read(vSpc,Ok,Msg)
    if(.not. Ok) &
    & call Stop_("SolModel_Read==> "//trim(Msg))

    call SolPhase_Read(vSpc,vSolModel,Ok,Msg)
    if(.not. Ok) &
    & call Stop_("SolPhase_Read==> "//trim(Msg))
    ! call SolPhase_Alloc ! allocate vSolFas

    ! call SolPhase_Init( & !
    ! & 1,          & !IN, model index
    ! & vSolModel,  & !IN, base of solution models
    ! & vSolFas(1), & !INOUT, solution phase
    ! & Ok,Msg)       !OUT
 
    ! call SolPhase_Model_Init(SolModel%Model,vSolModel,vSolFas(1),Ok,Msg)
    ! if(.not. Ok) call Stop_(trim(Msg))

  end if
  !-----------------------------------------------------/build vSolModel
  !
  call MixPhases_Alloc(vSpc,vMixModel,       vMixFas)
  !
  call Phases_Alloc_New(vSpc,vMixFas,vSolFas,   vFas)
  !
  call Condition_Read(iTP)
  !
  TdgK= Tref
  Pbar= Pref
  !
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
  call Species_Write_AquSize(vSpc)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Global_Build"
  !
end subroutine Global_Build_New

subroutine vSpecies_Rename(vSpc)
  use M_T_Species,only: T_Species,Species_Rename
  type(T_Species),intent(inout):: vSpc(:)
  !
  integer:: I
  !
  do I=1,size(vSpc)
    call Species_Rename(vSpc(I)%NamSp)
  end do
  !
end subroutine vSpecies_Rename

subroutine Condition_Read(iTP)
!--
!-- read the block CONDITIONS
!--
  use M_Trace,only: fHtm,iDebug,fTrc
  use M_Files,only: NamFInn,cTitle,File_Path
  use M_Files,only: DirOut,DirLog,DirDtbLog,DirDtbOut
  use M_Files,only: Files_Index_Close,Files_Index_Open
  use M_IoTools
  !
  integer,intent(out):: iTP
  !
  character(len=255):: L
  character(len=80) :: W
  logical:: EoL,Ok
  integer:: f,ios,N
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Condition_Read"
  !
  iTP= 0
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))

  Ok=.false.

  DoFile: do
    !
    read(f,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)

    case("ENDINPUT")
      exit  DoFile
    !
    case("CONDITIONS")
      Ok=.true.
      !
      LoopRead: do
        !
        read(f,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle LoopRead !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        select case(W)
        !
        case("ENDINPUT")            ; exit DoFile
        case("END","ENDCONDITIONS") ; exit LoopRead
        !
        case("TP.POINT","TPPOINT")
          !! call Stop_(trim(W)//" is Obsolete !! Give TdgC, Pbar in SYSTEM")
          call LinToWrd(L,W,EoL); call WrdToInt(W,iTP)
        case("TITLE")
          cTitle=trim(L)
        case("OUTPUT")
          call LinToWrd(L,W,EoL,"NO"); DirOut=File_Path(trim(W))
        case("DIRLOG")
          call LinToWrd(L,W,EoL,"NO"); DirLog=File_Path(trim(W))
        case("DTBDIROUT")
          call LinToWrd(L,W,Eol,"NO"); DirDtbOut=File_Path(trim(W))
        case("DTBDIRLOG")
          call LinToWrd(L,W,Eol,"NO"); DirDtbLog=File_Path(trim(W))
        case("DEBUG")
          call LinToWrd(L,W,EoL); call WrdToInt(W,N)
          if(N>0 .and. iDebug==0) iDebug=N
        !
        end select
        !
      end do LoopRead
    !
    end select
    !
  end do DoFile

  close(f)

  if(.not.Ok) then
    if(iDebug>0) write(fTrc,'(A)') "block CONDITIONS not found ???!!!"
    if(iDebug>2) print '(A)',"block CONDITIONS not found ???!!!"
  end if

  if(iDebug>0) then
    if(fHtm>0) call Files_Index_Close(fHtm)
    call Files_Index_Open(fHtm)
  end if

  if(iDebug>0) write(fTrc,'(A,/)') "</ Condition_Read"

end subroutine Condition_Read

end module M_Global_Build
