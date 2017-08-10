module M_DiscretModel_Test
  !.object describing "discretization" of a mixture phase
  !.to an array of pure phases
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  implicit none
  !
  private
  !
  public:: DiscretModel_Test
  !
contains

subroutine Stoikio_Test
  use M_IoTools,only: GetUnit
  use M_T_Species,only: T_Species
  use M_Global_Vars,only: vEle,vFas,vSpc,vMixModel
  use M_Global_Vars,only: vDiscretModel,vDiscretParam
  !
  integer:: f
  integer:: iSp,iEl,iDis
  type(T_Species):: S
  !
  call GetUnit(f)
  open(f,file="debug_discretmodel_stoik.log")
  !
  write(f,'(2(A,A1))',advance="no") "name",T_,"discri",T_
  write(f,'(3(A,A1))',advance="no") "disc%I",T_,"disc%iEl",T_,"disc%iDis",T_
  do iEl=1,size(vEle)
    write(f,'(A,A1)',advance="no") vEle(iEl)%NamEl,T_
  end do
  write(f,'(A)') "div"
  !
  do iSp=1,size(vSpc)
    !
    S= vSpc(iSp)
    iDis= S%iDiscret
    !
    if(iDis==0) cycle
    !
    write(f,'(A,A1)',advance="no") S%NamSp,T_
    write(f,'(A,A1)',advance="no") &
    & vDiscretModel(vDiscretParam(iDis)%iModel)%Name,T_
    !
    write(f,'(3(I3,A1))',advance="no") &
    & vDiscretParam(iDis)%i,T_, &
    & vDiscretParam(iDis)%j,T_, &
    & vDiscretParam(iDis)%k,T_
    !
    do iEl=1,size(vEle)
      write(f,'(I3,A1)',advance="no")  S%vStoikio(iEl),T_
    end do
    write(f,'(I3,A1)') S%vStoikio(0)
    !
  end do
  !
  close(f)
  !
  if(iDebug>1) then
    print *,">> results in debug_discretmodel_stoik.log !!!"
    call Pause_
  end if
end subroutine Stoikio_Test

subroutine DiscretModel_Test(TdgK,Pbar)
  use M_Global_Vars,only: vEle,vFas,vSpc,vMixModel
  use M_Global_Vars,only: vDiscretModel,vDiscretParam
  use M_DiscretModel_Tools
  !
  real(dp),intent(in):: TdgK,Pbar
  !
  integer:: i
  !
  if(idebug>1) write(fTrc,'(/,A)') "< DiscretModel_Test"
  !
  call DiscretModel_Test_Init(TdgK,Pbar)
  !
  if(iDebug>1) then
    !
    print *,"<============================================== SPECIES =="
    do i=1,size(vSpc)
      print *,vSpc(i)%NamSp
    end do
    call Pause_
    !
    print *,"<=========================================== MIX.MODELS =="
    do i=1,size(vMixModel)
      print *,vMixModel(i)%Name
    end do
    call Pause_
    !
    print *,"<======================================= DISCRET.MODELS =="
    do i=1,size(vDiscretModel)
      print *,vDiscretModel(i)%Name
    end do
    call Pause_
    !
    print *,"<=============================================== PHASES =="
    do i=1,size(vFas)
      print *,vFas(i)%NamFs
    end do
    call Pause_
    !
    print *,"<======================================= DISCRET.PARAMS =="
    do i=1,size(vDiscretParam)
      print *, &
      & vDiscretModel(vDiscretParam(i)%iModel)%Name, &
      & vDiscretParam(i)%I, &
      & vDiscretParam(i)%J, &
      & vDiscretParam(i)%K
    end do
    call Pause_
    !
  end if
  !
  call DiscretSpecies_Stoikio_Calc( &
  & vEle,          & !IN
  & vMixModel,     & !IN
  & vDiscretModel, & !IN
  & vDiscretParam, & !IN
  & vSpc)            !INOUT with updated stoichio's for discretized species
  !
  call Stoikio_Test
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ DiscretModel_Test"
  !
end subroutine DiscretModel_Test

subroutine DiscretModel_Test_Init(TdgK,Pbar)
  use M_T_Species,   only: T_Species,Species_Stoikio_Calc,Species_Append
  use M_Dtb_Const,   only: T_CK,Tref,Pref
  use M_Dtb_Read,    only: Dtb_Read_Species
  use M_Dtb_Calc,    only: SpeciesDtb_ToSpecies
  use M_Element_Read,only: Element_Read_Redox,Element_Read_Entropy
  use M_Global_Tools,only: Global_TP_Update
  use M_DiscretModel_Read
  use M_DiscretModel_Tools
  use M_Global_Alloc
  !
  use M_Global_Vars, only: vEle,vSpc,vSpcDtb,vFas
  use M_Global_Vars, only: vMixModel,vMixFas
  use M_Global_Vars, only: vDiscretModel,vDiscretParam
  use M_Dtb_Vars,    only: DtbLogK_vTPCond
  !
  real(dp),intent(in):: TdgK,Pbar
  !
  integer:: N
  logical:: fOk
  type(T_Species),allocatable:: vSpcTmp(:)
  !
  call Elements_Alloc_forDtb(vEle,N) !---------------------> builds vEle
  !
  call Element_Read_Redox(vEle)
  call Element_Read_Entropy(vEle)
  !
  call Dtb_Read_Species(vEle) !--------> read databases -> build vDtbXxx
  !
  call SpeciesDtb_Alloc(vSpcDtb,N) !---------> from vDtb*, build vSpcDtb
  !
  if(N<1) return
  !
  deallocate(vSpc)
  allocate(vSpc(1:N))
  !
  !! call Species_Alloc_forDtb(vEle) !-> read databases -> vSpc,tFormula
  call SpeciesDtb_ToSpecies(vSpcDtb,vSpc)
  !
  call Species_Stoikio_Calc(vEle,vSpc,fOk)
  !
  call MixModels_Alloc(vSpc,vMixModel) !--> build MixModel dtb vMixModel
  !
  call DiscretModel_Read(vMixModel) !-----------> allocate vDiscretModel
  !
  N= size(vDiscretModel)
  if(N>0) then
    call DiscretParam_Alloc(vDiscretModel,   vDiscretParam)
    !
    allocate(vSpcTmp(size(vDiscretParam)))
    !
    call DiscretParam_Init( &
    & vEle,vSpc,vMixModel,vDiscretModel, &
    & vDiscretParam,vSpcTmp) !-> build vSpcDiscret
    !
    call Species_Append(vSpcTmp,     vSpc) !-> new vSpc !!!
    !
    deallocate(vSpcTmp)
    !
  end if
  !
  !-------------------------------read phase compositions, build vMixFas
  call MixPhases_Alloc(vSpc,vMixModel,  vMixFas) 
  !
  call Phases_Alloc(vSpc,vMixFas, vFas) !-------------------> build vFas
  !
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
end subroutine DiscretModel_Test_Init

end module M_DiscretModel_Test
