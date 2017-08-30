module M_GEM_Theriak
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  use M_Kinds
  use M_Trace,       only: iDebug,fTrc,T_,Stop_,Pause_
  use M_T_Phase,     only: T_Phase
  use M_T_MixPhase,  only: T_MixPhase
  use M_T_MixModel,  only: T_MixModel,MaxPole,Mix_Molecular
  use M_GEM_Vars,    only: T_SavPhase
  !
  implicit none
  !
  private
  !
  public:: GEM_Theriak
  public:: GEM_Theriak_Path
  !
  type:: T_GemPhase
    character(len=23):: NamFs
    integer :: nC
    real(dp):: vStoik(1:12) ! vStoik(1:MaxEle)
    integer :: iModel
    real(dp):: vXPole(1:MaxPole)
    real(dp):: Grt
    real(dp):: Mole
  end type T_GemPhase
  !
  ! type:: T_Phase
  !   character(len=23):: NamFs
  !   integer :: iSpc != index of species in vSpc, in case of pure phase 
  !   integer :: iMix != index of mixture in vMixFas, in case of mixture phase
  !   integer :: iSol != index of solution in vSolFas, in case of solution phase
  !   real(dp):: Grt,VolM3,WeitKg
  !   real(dp):: Mole
  ! end type T_Phase
  !
  !----------------------------------------------------private variables
  integer:: F1,F2,FF
  !
  type(T_SavPhase),allocatable:: vSavPhase(:)
  type(T_SavPhase),allocatable:: tResultMix(:,:)
  !
  real(dp),allocatable:: tResult(:,:)
  real(dp),allocatable:: vFasMolPur(:)
  real(dp),allocatable:: vFasMolMix(:)
  logical, allocatable:: vFasIsPresent(:)
  logical, allocatable:: vModelConvex(:)
  !--------------------------------------------------------------------/
  !
  !---------------------------------------------------private parameters
  real(dp),parameter:: MixMinim_TolX= 1.D-4
  real(dp),parameter:: GEM_G_Iota= 1.0D-9 !1.0D-9
  integer, parameter:: GEM_IterMax= 500
  !
  logical:: WarmRestart= .false.  !!.true.   !!
  logical:: WriteTrace=  .true.   !!.false.  !!
  !--------------------------------------------------------------------/
  !
contains

subroutine GEM_Theriak
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  use M_IoTools,     only: GetUnit
  use M_Files,       only: DirOut
  use M_Simplex_Vars,only: tSimplex
  use M_Simplex_Vars,only: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  use M_Global_Vars, only: vFas
  use M_GEM_Vars,    only: vCpnGEM,tStoikioGEM
  !
  !--- for Global_TP_Update
  use M_Global_Tools,only: Global_TP_Update
  use M_Global_Vars, only: vSpcDtb,vSpc,vMixModel
  use M_Global_Vars, only: vDiscretModel,vDiscretParam,vMixFas
  use M_GEM_Vars,    only: TdgK,Pbar
  !---/
  !
  type(T_SavPhase):: SavPhaseZero
  integer:: iError
  integer:: I,nC,nFpur,nMix
  integer :: fo !!JM!!216-10
  real(dp):: x
  !---------------------------------------------------------------------
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< GEM_Theriak"
  !
  F1= 0
  F2= 0
  FF= 0
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  nC= size(vCpnGEM)
  nFpur= size(vFas)
  nMix= size(vMixModel)
  !
  !----------------------------------------------------------allocations
  allocate(vFasMolPur(1:nFpur))          ;  vFasMolPur(:)= Zero
  allocate(vFasIsPresent(1:nFpur+nMix))  ;  vFasIsPresent=.false.
  !
  if(iDebug>2) then
    print *,"nMix=",nMix
    call pause_
  end if
  !
  allocate(vFasMolMix(nMix))
  allocate(vSavPhase(nMix))
  allocate(vModelConvex(nMix))
  !
  if(nMix>0) then
    !
    vFasMolMix(:)= Zero
    vSavPhase(:)= SavPhaseZero
    !
    vModelConvex(:)=.false.
    !
  end if
  !---------------------------------------------------------/allocations
  
  !-----------------------------------------------------------open files
  if(iDebug>2) then
    if(nMix>0) then
      call Write_Log_Entete(vFas)
      if(iDebug>3) then
        call GetUnit(ff)
        open(ff,file='log_theriak_minimize.log')
      end if
    end if
  end if
  !----------------------------------------------------------/open files
  
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  
  !--------------------------------------------------------------SIMPLEX
  call Simplex_GEM(iError)
  !-------------------------------------------------------------/SIMPLEX
  !
  !if(iError/=0) then
  !  call Error_Show(iError)
  !else
  !  do I=1,nFpur
  !    if(vFasMolPur(I)>Zero) write(11,'(G15.6,2X,A)') &
  !    & vFasMolPur(I),vFas(I)%NamFs
  !  end do
  !  if(nMix>0) then
  !    do I=1,nMix
  !      if(vFasMolMix(I)>Zero) write(11,'(G15.6,2X,A)') &
  !      & vFasMolMix(I),vMixModel(I)%Name
  !    end do
  !  end if
  !end if

  !modif-----------------------------------------------------JM--2016-10
  call GetUnit(fo)
  ! open(fo,file=trim(DirOut)//"_gem.out")
  open(fo,file="tmp_gem.tab")
  !
  if(iError/=0) then
    write(fo,'(A)') "error"
    call Error_Show(iError)
  else
    write(fo,'(4(A,A1))') &
    & "PHASE",        T_, &
    & "MOLES",        T_, &
    & "VOLUME_M3",    T_, &
    & "DENSITY_KG_M3",T_
    !
    do I=1,nFpur
      if(vFasMolPur(I)>Zero) write(fo,'(A,A1,3(G15.6,A1))') &
      & trim(vFas(I)%NamFs),         T_, &
      & vFasMolPur(I),               T_, &
      & vFasMolPur(I)*vFas(I)%VolM3, T_, &
      & vFas(I)%WeitKg/vFas(I)%VolM3,T_
    end do
    !
    x= 0._dp
    do I=1,nMix
      if(vFasMolMix(I)>Zero) write(fo,'(A,A1,3(G15.6,A1))') &
      & trim(vMixModel(I)%Name),T_, &
      & vFasMolMix(I),          T_, &
      & x,                      T_, &
      & x,                      T_
    end do
  end if
  close(fo)
  !---------------------------------------------------------/JM--2016-10

  if(F1>0) close(F1)
  if(F2>0) close(F2)
  if(FF>0) close(ff)

  deallocate(vFasMolPur)
  deallocate(vFasIsPresent)
  !
  deallocate(vFasMolMix)
  deallocate(vSavPhase)
  deallocate(vModelConvex)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ GEM_Theriak"

  return
end subroutine GEM_Theriak

subroutine GEM_Theriak_Path
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!-- changing step by step the amount of components' mole numbers
!--
  use M_IoTools,     only: GetUnit
  use M_Dtb_Const,   only: T_CK
  use M_Path_Read,   only: Path_ReadMode, Path_ReadParam_new
  use M_Simplex_Vars,only: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  use M_Files,       only: NamFInn
  use M_Files,       only: DirOut
  use M_IoTools,     only: GetUnit
  !
  use M_Global_Vars, only: vFas
  use M_GEM_Vars,    only: vCpnGEM,tStoikioGEM
  use M_GEM_Vars,    only: TdgK,Pbar
  !
  use M_Path_Vars,   only: vTPpath,vLPath,tPathData,DimPath
  use M_Path_Vars,   only: Path_Vars_Clean
  use M_Simplex_Vars,only: tSimplex
  use M_GEM_Write,   only: GEM_Write_Phases,GEM_Write_Mixtures
  !
  !--for Global_TP_Update
  use M_Global_Tools,only: Global_TP_Update
  use M_Global_Vars, only: vSpcDtb,vSpc,vMixModel
  use M_Global_Vars, only: vDiscretModel,vDiscretParam,vMixFas
  !---------------------------------------------------------------------
  !
  type(T_SavPhase):: SavPhaseZero
  integer :: iPath,I,J
  integer :: iError
  integer :: nC,nFpur,nMix
  real(dp):: TdgK0,Pbar0
  !
  logical :: Ok
  character(len=3) :: PathMod3
  character(len=80):: Msg
  !
  logical, allocatable:: vSimplex_Ok(:)
  !---------------------------------------------------------------------
  
  if(iDebug>0) write(fTrc,'(/,A)') "< GEM_Theriak_Path"
  !
  nC= size(vCpnGEM)
  nFpur= size(vFas) ! + size(vMixModel)
  nMix= size(vMixModel)
  !
  F1= 0
  F2= 0
  FF= 0
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  !---------------------------------read path parameters from PATH block
  call Path_ReadMode(NamFInn,PathMod3,Ok,Msg)  !out
  !
  if(PathMod3 /= "CHG" .and. &
  &  PathMod3 /= "GRD" ) then
    Msg= "Global equilibrium paths: only in CHANGE or GRID modes"
    Ok= .false.
    print *,trim(Msg)
    return
  end if
  !
  if(PathMod3 == "CHG") call Path_ReadParam_new( &
  & NamFInn,  &
  & PathMod3, &
  & vCpnGEM, &
  & TdgK,Pbar, &
  & Ok,Msg)
  !
  if(.not. Ok) then
    print *,trim(Msg)
    return
  end if
  !--------------------------------/read path parameters from PATH block
  
  !----------------------------------------------------------allocations
  allocate(vFasMolPur(1:nFpur))          ;  vFasMolPur(:)= Zero
  allocate(vFasIsPresent(1:nFpur+nMix))  ;  vFasIsPresent=.false.
  !
  allocate(vFasMolMix(nMix))
  allocate(vSavPhase(nMix))
  allocate(vModelConvex(nMix))
  !
  if(allocated(tResult)) deallocate(tResult)
  allocate(tResult(1:nC+nFpur+nMix+2,1:dimPath))  ; tResult=Zero
  if(nMix>0) allocate(tResultMix(1:nMix,1:dimPath))
  !
  allocate(vSimplex_Ok(1:dimPath))  ;  vSimplex_Ok=.false.
  !
  if(nMix>0) then
    !
    vFasMolMix(:)= Zero
    vSavPhase(:)= SavPhaseZero
    !
    vModelConvex(:)=.false.
    !
  end if
  !---------------------------------------------------------/allocations
  !
  !-----------------------------------------------------------open files
  if(nMix>0) then
    if(iDebug>2) call Write_Log_Entete(vFas)
    !
    if(iDebug>3) then
      call GetUnit(ff)
      open(ff,file='log_theriak_minimize.log')
    end if
  end if
  !----------------------------------------------------------/open files
  !
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  !
  !do i=1,size(vFas)
  !  print *, vFas(i)%VolM3*1.D6," = ",trim(vFas(i)%NamFs)
  !end do
  !pause

  do iPath=1,dimPath
    !
    TdgK0= TdgK
    Pbar0= Pbar
    !
    !----------------------if change T,P update T,P-dependent properties
    TdgK= vTPpath(iPath)%TdgC +T_CK
    Pbar= vTPpath(iPath)%Pbar
    !
    if(TdgK0 /= TdgK .or. Pbar0 /= Pbar) &
    & call Global_TP_Update( &
    & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
    & vSpc,vMixModel,vMixFas,vFas)
    !------------------------------------------------------------------/
    !
    !-------------------------------------------------system composition
    ! print *,"system composition:"
    do J=1,nC
      if(vLPath(J)) then
        vCpnGEM(J)%Mole= tPathData(J,iPath)
      end if
    enddo
    ! call pause_
    !------------------------------------------------------------------/
    !
    tResult(1,     iPath)= TdgK -T_CK
    tResult(2,     iPath)= Pbar
    tResult(3:nC+2,iPath)= vCpnGEM(1:nC)%Mole
    !
    !--------------------------------------------------------SIMPLEX_GEM
    call Simplex_GEM(iError)
    !-------------------------------------------------------/SIMPLEX_GEM
    !
    !------------------------store results in tables tResult, tResultMix
    if(iError/=0) then
      !
      vSimplex_Ok(iPath)= .false.
      call Error_Show(iError)
      !
    else
      !
      vSimplex_Ok(iPath)= .true.
      !
      tResult(nC+3:nC+nFpur+2,iPath)= vFasMolPur(1:nFpur)
      if(nMix>0) then
        tResult(nC+nFpur+3:nC+nFpur+nMix+2,iPath)= vFasMolMix(1:nMix)
        tResultMix(1:nMix,iPath)= vSavPhase(1:nMix)
        do i=1,nMix
          ! print *,vSavPhase(i)%vVol0(1)
          if(vSavPhase(i)%nFas > 1) call Check_vXMean(i)
        enddo
        ! call pause_ !!
      end if
      !
      ! call Path_StoreResult(iPath,vCpnGEM(:)%Mole,TdgK,Pbar)
      !
    end if
    !-----------------------/store results in tables tResult, tResultMix
    !
    if(iDebug==1) print *,iPath
    if(iDebug>1) then
      print *, &
      & "============================================================", &
      & iPath
      if(iDebug>2) call pause_
    end if
    !
  enddo
  !
  !--------------------------------------------------write result tables
  call GEM_Write_Phases( &
  & DimPath, &
  & vFasIsPresent, &
  & vSimplex_Ok, &
  & tResult, &
  & tResultMix)
  !
  if(size(vMixModel)>0) call GEM_Write_Mixtures( &
  & DimPath, &
  & vFasIsPresent, &
  & vSimplex_Ok,&
  & MixMinim_TolX, &
  & tResult, &
  & tResultMix)
  !-------------------------------------------------/write result tables
  !
  deallocate(vSimplex_Ok)
  deallocate(tResult)
  deallocate(vFasMolPur)
  deallocate(vFasIsPresent)
  deallocate(vFasMolMix)
  deallocate(vSavPhase)
  deallocate(vModelConvex)
  !
  if(nMix>0) deallocate(tResultMix)
  !
  call Path_Vars_Clean
  !
  if(F1>0) close(F1)
  if(F2>0) close(F2)
  if(ff>0) close(ff)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ GEM_Theriak_Path"
  !
contains

subroutine Check_vXMean(i)
!--
!-- for a mixing model with more than one phase in stable assemblage
!-- compute the Gibbs for the mean composition of the different phases.
!-- if the Gibbs is not close to zero, then there is unmixing
!--
  integer,intent(in):: i

  integer :: J,nP,ff
  real(dp):: G
  real(dp),allocatable:: vX(:)

  J=  vSavPhase(i)%iModel
  nP= vMixModel(J)%NPole
  allocate(vX(1:nP))
  vX(:)= Zero
  !
  do ff=1,vSavPhase(i)%nFas
    vX(:)= vX(1:nP) &
    &    + vSavPhase(i)%vMole(ff) *vSavPhase(i)%tXPole(ff,1:nP)
  enddo
  vX(:)= vX(:) /SUM(vX(:))
  !
  G= MixPhase_Grt( &
  & TdgK,Pbar,    &
  & nP, &
  & vMixModel(J), &
  & vSavPhase(i)%vGrt0(1:nP), &
  & vX(1:nP))
  !
  deallocate(vX)

  !~ print *,"G(vXMean)..",vMixModel(J)%Name, G  ;  pause

  return
end subroutine Check_vXMean

end subroutine GEM_Theriak_Path

subroutine Simplex_GEM(iError)
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  use M_Simplex_Calc,only: Simplex_Calc
  use M_Simplex_Vars,only: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  use M_Global_Vars, only: vFas,vMixModel
  use M_GEM_Vars,    only: vCpnGEM,tStoikioGEM
  use M_GEM_Vars,    only: TdgK,Pbar
  !
  use M_Simplex_Vars,only: iPosV,tSimplex
  !
  integer, intent(out):: iError
  !
  real(dp),        allocatable:: tStoikio(:,:)
  type(T_Phase),   allocatable:: vFasPur(:)
  type(T_GemPhase),allocatable:: vFas0(:)
  type(T_MixPhase),allocatable:: vMixFas0(:)
  type(T_SavPhase),allocatable:: vMixFas_Xpole_Init(:)
  type(T_SavPhase):: SavPhaseZero
  !
  ! type(T_MixModel):: MM
  !
  integer :: nC,nF
  integer :: nFpur,nFmix,nMix
  integer :: nFmixSpl
  integer :: Iter
  integer :: I,K
  real(dp):: G_Iota
  ! character(len=30):: sFMT
  !
  G_Iota= GEM_G_Iota
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  nC=    size(vCpnGEM)
  nFpur= size(vFas) !== at this point, only pure phases are in vFas ==
  nMix=  size(vMixModel)
  !
  nFmix= nMix
  !
  !-----------------------------------------------------------Allocation
  allocate(tStoikio(1:nFpur+nFmix+nC,1:nC))
  allocate(vFasPur(1:nFpur))
  allocate(vFas0(1:nFpur+nFmix+nC))
  !-- vMixFas0: One mixture phase for Each mixture model
  allocate(vMixFas0(nFmix))
  allocate(vMixFas_Xpole_Init(nFmix))
  !----------------------------------------------------------/Allocation
  !
  tStoikio(1:nFpur,:)= tStoikioGEM(1:nFpur,:)
  vFasPur(1:nFpur)= vFas(1:nFpur)
  !
  do I=1,nFpur
    vFas0(I)%NamFs=     vFas(I)%NamFs
    vFas0(I)%nC=        nC
    vFas0(I)%vStoik(1:nC)= tStoikioGEM(I,1:nC)
    vFas0(I)%iModel=    0 != no mixing model, because phase is pure
    vFas0(I)%Grt=       vFas(I)%Grt
    vFas0(I)%Mole=      Zero
  enddo
  !
  if(nMix>0) then
    do K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    enddo
    call MixPhase_XPole_Init(vMixModel,vMixFas_Xpole_Init)
  end if

  !------------------------------------------------------initial simplex
  nFmixSpl= 0
  nF= nFpur
  !
  call Simplex_Vars_Alloc(nC,nF)
  !
  tSimplex(1:nC, 0)=  vCpnGEM(1:nC)%Mole
  tSimplex(0,    1:nF)= -vFas0(1:nF)%Grt
  tSimplex(1:nC, 1:nF)= -TRANSPOSE(tStoikio(1:nF,1:nC))
  !
  call Simplex_Calc(iError)
  !
  if(iError/=0) then
    deallocate(tStoikio)
    deallocate(vFasPur)
    deallocate(vFas0)
    ! if(nMix>0) deallocate(vMixFas0)
    ! if(nMix>0) deallocate(vMixFas_Xpole_Init)
    deallocate(vMixFas0)
    deallocate(vMixFas_Xpole_Init)
    return
  end if
  !
  !-----------------compute new free energies of formation of all phases
  !---------------------------------------from current stable assemblage
  call Gibbs_Change( &
  & nC,nF,IPOSV, &
  & vFas0)
  !--------------------------------------------------------------------/
  !
  ! if(iDebug>2) &
  ! & call GEM_ShowResults(nC,vMixModel,vFas0,IPOSV,tSimplex(:,0))
  ! if(iDebug>2) pause
  !
  !-----------------------------------------------------/initial simplex

  if(nMix==0) then
    !
    if(iDebug>2) call ShowResults
    call StoreResults
    !
  else
    !
    Iter= 0
    do
      !
      Iter= Iter+1
      !
      !-------------------------------calc' minimals G's of all mixtures
      vFasPur(1:nFpur)%Grt= vFas0(1:nFpur)%Grt
      !
      call Mixture_Minimize( &
      & TdgK,Pbar, &
      & vFasPur,vMixModel,  &
      & vMixFas_Xpole_Init, &
      & vMixFas0)
      !----------------------------------------------------------------/
      if(ALL(vMixFas0(:)%Grt >=-G_Iota)) then
        if(iDebug>1) call ShowResults
        call StoreResults
        exit
      end if
      !
      !-------- to the current table of G and stoichiometry for simplex,
      !----------------------- append all mixture phases with negative G
      !---------------------------------- as phases of fixed composition
      !------------- (with the composition computed in Mixture_Minimize)
      nF= nFpur + nFmixSpl
      call GEM_AddMixtures( &
      & nC, G_Iota, vMixModel, vMixFas0, &
      & vFas0, nF)
      !
      do I=nFpur+1,nF
        tStoikio(I,1:nC)= vFas0(I)%vStoik(1:nC)
      enddo
      !
      call Simplex_Vars_Clean
      call Simplex_Vars_Alloc(nC,nF)
      !
      tSimplex(1:nC,0   )=  vCpnGEM(1:nC)%Mole
      tSimplex(0,   1:nF)= -vFas0(1:nF)%Grt
      tSimplex(1:nC,1:nF)= -TRANSPOSE(tStoikio(1:nF,1:nC))
      !
      call Simplex_Calc(iError)
      !
      if(iError/=0) then
        deallocate(tStoikio)
        deallocate(vFasPur)
        deallocate(vFas0)
        ! if(nMix>0) deallocate(vMixFas0)
        ! if(nMix>0) deallocate(vMixFas_Xpole_Init)
        deallocate(vMixFas0)
        deallocate(vMixFas_Xpole_Init)
        return
      end if
      !
      if(iDebug>1) call ShowResults
      call StoreResults
      !
      !-------------compute new free energies of formation of all phases
      !-----------------------------------from current stable assemblage
      call Gibbs_Change( &
      & nC,nF,IPOSV, &
      & vFas0)
      !--/
      !
      if(iDebug>2) then
        write(92,'(/,A,I3)') "vFas0%Grt, iter= ", Iter
        do I=1,nF
          write(92,'(4X,I3,2X,G15.6,2X,A)') &
          & I,vFas0(I)%Grt,trim(vFas0(I)%NamFs)
        enddo
      end if
      !
      !--------------------------------to vFas(1:nFpur) (= pure phases),
      !---append mixture phases that belong to current stable assemblage
      !-------------------------------------found from last simplex call
      call GEM_AppendSPLMixtures(nC,nFpur,iPosV,vFas0,nFmixSpl)
      !
      !call Simplex_Vars_Clean
      !
      if(Iter > GEM_IterMax) exit
      !
    enddo
    !
    if(iDebug>2) print '(A,I3)',"Iter= ",Iter
    !
  end if
  !
  call GEM_SaveResults(nFpur,Iter,vMixModel)
  !
  call Simplex_Vars_Clean
  !--------------------------------------------------------DesAllocation
  deallocate(tStoikio)
  deallocate(vFasPur)
  deallocate(vFas0)
  ! if(nMix>0) deallocate(vMixFas0)
  ! if(nMix>0) deallocate(vMixFas_Xpole_Init)
  deallocate(vMixFas0)
  deallocate(vMixFas_Xpole_Init)
  !-------------------------------------------------------/DesAllocation
  !
contains

  subroutine ShowResults
    integer:: I,J,K,P,M
    
    write(6,'(A)') "pure="
    do K=1,nC
      I= IPOSV(K)
      if(vFas0(I)%iModel==0) &
      & write(6,'(4X,G15.6,2X,A)') tSimplex(K,0), trim(vFas0(I)%NamFs)
    enddo
    
    write(6,'(A)') "MIX="
    do K=1,nC
      I= IPOSV(K)
      if(vFas0(I)%iModel/=0) then
        write(6,'(4X,2I3,2X,G15.6,A)',advance="NO") &
        & I,K,tSimplex(K,0)," X="
        J= vFas0(I)%iModel
        do P=1,vMixModel(J)%NPole
          write(6,'(F7.3,1X)',advance="NO") vFas0(I)%vXPole(P)
        enddo
        write(6,'(2X,A)') trim(vFas0(I)%NamFs)
      end if
    enddo
    write(6,*)
    
    return
  end subroutine ShowResults

  subroutine StoreResults
    integer:: I,J,K,P
    integer:: nP

    vFasMolPur(:)= Zero
    !
    vSavPhase(:)= SavPhaseZero
    do I=1,size(vSavPhase)
      vSavPhase(I)%iModel= I
    enddo

    do K=1,nC
      !
      I= IPOSV(K)
      J= vFas0(I)%iModel
      !
      if(vFas0(I)%iModel==0) then
        !
        vFasMolPur(I)= tSimplex(K,0)
        !
      else
        !
        vSavPhase(J)%iModel= J
        vSavPhase(J)%nFas= vSavPhase(J)%nFas +1
        vSavPhase(J)%vMole(vSavPhase(J)%nFas)= tSimplex(K,0)
        nP= vMixModel(J)%NPole
        do P=1,nP
          vSavPhase(J)%tXPole(vSavPhase(J)%nFas,P)= vFas0(I)%vXPole(P)
          vSavPhase(J)%vGrt0(P)= vFasPur(vMixModel(J)%vIPole(P))%Grt
          vSavPhase(J)%vVol0(P)= vFasPur(vMixModel(J)%vIPole(P))%VolM3
          !print *,"StoreResults",vSavPhase(J)%vVol0(P)
        enddo
        !
      end if
      !
    enddo

    !--------------------------------------------------------------trace
    if(iDebug>2) then
      do K=1,nC
        I= IPOSV(K)
        if(vFas0(I)%iModel/=0) then
          write(93,'(A,A1,G15.6,A1)',advance="NO") &
          & trim(vFas0(I)%NamFs),T_, &
          & tSimplex(K,0),       T_
          J= vFas0(I)%iModel
          do P=1,vMixModel(J)%NPole
            write(93,'(G15.6,A1)',advance="NO") vFas0(I)%vXPole(P),T_
          enddo
        end if
      enddo
      write(93,*)
    end if
    !-------------------------------------------------------------/trace

    return
  end subroutine StoreResults

end subroutine Simplex_GEM

real(dp) function MixPhase_Grt(TdgK,Pbar,nP,MM,vMu0rt,vX)
  use M_T_MixModel,only: MixModel_Activities
  !---------------------------------------------------------------------
  real(dp),        intent(in):: TdgK,Pbar
  integer,         intent(in):: nP
  type(T_MixModel),intent(in):: MM        ! mixing model
  real(dp),        intent(in):: vMu0rt(:) !
  real(dp),        intent(in):: vX(:)     ! phase composition
  !---------------------------------------------------------------------
  real(dp):: vLGam(nP),vLIdeal(nP),vLnAct(nP)
  logical :: vLPole(nP)
  real(dp):: G
  integer :: i
  logical :: Ok
  character(len=30):: Msg
  !---------------------------------------------------------------------
  vLPole(:)= (vX(:)>Zero) ! .and. MM%vHasPole(:)

  call MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MM,        & ! in
  & vX,        & ! in
  & vLPole,    & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & !
  & vLIdeal,   & !
  & vLnAct)      !
  !
  G= Zero
  do i=1,nP
    if(vLPole(i)) &
    ! vMu0rt(i)= vFasPur(MM%vIPole(i))%Grt
    & G= G &
    &  + vX(i) *(vMu0rt(i) + vLnAct(i))
  end do
  !
  MixPhase_Grt= G

  return
end function MixPhase_Grt

subroutine GEM_SaveResults(nFpur,Iter,vMixModel)
  integer,         intent(in):: nFpur
  integer,         intent(in):: Iter
  type(T_MixModel),intent(in):: vMixModel(:)
  !
  integer :: I,J
  
  !--- pure phases --
  do I=1,nFpur
    if(vFasMolPur(I)>Zero) vFasIsPresent(I)= .true.
  enddo
  
  !--- mixtures --
  vFasMolMix(:)= Zero
  do I=1,size(vSavPhase)
    if(vSavPhase(I)%nFas >0) then
      vFasIsPresent(nFpur +I)= .true.
      do J=1,vSavPhase(I)%nFas
        vFasMolMix(I)= vFasMolMix(I) +vSavPhase(I)%vMole(J)
      enddo
    end if
  enddo
  
  return
end subroutine GEM_SaveResults

subroutine GEM_AppendSPLMixtures(nC,nFpur,iPosV,vFas,nF)
!--
!-- to vFas(1:nFpur) (- pure phases),
!-- append mixture phases that belong to current stable assemblage
!-- found from last simplex call
!--
  integer,         intent(in)   :: nC, nFpur
  integer,         intent(in)   :: IPOSV(:)
  type(T_GemPhase),intent(inout):: vFas(:)
  integer,         intent(inout):: nF
  !
  type(T_GemPhase),allocatable:: vFas1(:)
  integer:: I,K
  !
  allocate(vFas1(size(vFas)))
  vFas1= vFas
  nF= 0
  do K=1,nC
    I= IPOSV(K)
    if(vFas1(I)%iModel/=0) then
      nF= nF+1
      vFas(nFpur+nF)= vFas1(I)
    end if
  enddo
  deallocate(vFas1)
  !
end subroutine GEM_AppendSPLMixtures

subroutine GEM_AddMixtures( &
& nC, G_Iota, vMixModel, vMixFas, &
& vFas, nFas)
!--
!-- to the current table of G and stoichiometry for simplex calculations,
!-- append all mixture phases with negative G
!-- as phases of fixed composition
!-- (with the composition that makes G minimal, computed in Mixture_Minimize)
!--
  integer,         intent(in)   :: nC
  real(dp),        intent(in)   :: G_Iota
  type(T_MixModel),intent(in)   :: vMixModel(:)
  type(T_MixPhase),intent(in)   :: vMixFas(:)
  type(T_GemPhase),intent(inout):: vFas(:)
  integer,         intent(inout):: nFas
  !
  type(T_MixModel):: MM
  integer :: iMF,iMdl,iMix,iCp,N

  iMF= nFas
  !
  do iMix=1,size(vMixFas)
  
    if(vMixFas(iMix)%Grt <-G_Iota) then
      !
      iMdl= vMixFas(iMix)%iModel
      iMF=  iMF +1
      MM=   vMixModel(iMdl)
      N=    MM%NPole
      !
      vFas(iMF)%NamFs=       vMixFas(iMix)%Name
      vFas(iMF)%vXPole(1:N)= vMixFas(iMix)%vXPole(1:N)
      vFas(iMF)%iModel=      iMdl
      vFas(iMF)%Grt=         vMixFas(iMix)%Grt
      !
      do iCp=1,nC ! stoichiometry of mixture vs components
        vFas(iMF)%vStoik(iCp)= &
        & SUM( vFas(MM%vIPole(1:N))%vStoik(iCp) &
        &    * vMixFas(iMix)%vXPole(1:N) )
      enddo
      !
      if(iDebug>2) print *,"ADDED= ",trim(vFas(iMF)%NamFs)
      !
    end if
  
  enddo
  !
  nFas= iMF
  !
end subroutine GEM_AddMixtures

subroutine Error_Show(iErr)
  integer,intent(in):: iErr
  select case(iErr)
    case(1)
      write(*,'(A)') "Unbounded Objective Function"
    case(-1)
      write(*,'(A)') "No Solutions Satisfy Constraints Given"
    case(2)
      write(*,'(A)') "Max Iteration Reached !!!"
  end select
end subroutine Error_Show

subroutine Gibbs_Change( &
& nC,nF,IPOSV, &
& vFas0)
!--
!-- use stable phase assemblage as the new component set
!-- and recalculate free energies of formation of all species,
!-- as free enegies of formation from these components
!--
  use M_Numeric_Mat,only: LU_Decomp, LU_BakSub
  !---------------------------------------------------------------------
  integer, intent(in):: nC,nF
  integer, intent(in):: IPOSV(:)
  type(T_GemPhase),intent(inout):: vFas0(:)
  !---------------------------------------------------------------------
  real(dp),allocatable:: tTransform(:,:)
  integer, allocatable:: vIndx(:)
  real(dp),allocatable:: vY(:)
  real(dp),allocatable:: vGrt(:),vGrt0(:)
  !
  integer :: I,K
  logical :: bSingul
  real(dp):: D
  !---------------------------------------------------------------------
  allocate(tTransform(1:nC,1:nC))
  allocate(vIndx(1:nC),vY(1:nC))
  allocate(vGrt(1:nF))
  allocate(vGrt0(1:nC))
  !
  !---------------------------------build stoikio table of stable phases
  do K=1,nC
    I= IPOSV(K) ! index of stable phase in phase list
    tTransform(1:nC,K)= vFas0(I)%vStoik(1:nC)
    vGrt0(K)= vFas0(I)%Grt !-> free energy of stable phase
  enddo
  !--------------------------------------------------------------------/
  !
  call LU_Decomp(tTransform,vIndx,D,bSingul)
  if(bSingul) call Stop_("SINGUL IN Gibbs_Change")
  !
  do I=1,nF
    vY(1:nC)= vFas0(I)%vStoik(1:nC)
    !-----------------------Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
    call LU_BakSub(tTransform,vIndx,vY)
    vGrt(I)= vFas0(I)%Grt - SUM(vY(:)*vGrt0(:))
  enddo
  vFas0(1:nF)%Grt= vGrt(1:nF)
  !
  deallocate(tTransform)
  deallocate(vGrt,vGrt0)
  deallocate(vIndx,vY)
  !
end subroutine Gibbs_Change

subroutine MixPhase_XPole_Init(vMixModel,vSavFas)
  type(T_MixModel),intent(in):: vMixModel(:)
  type(T_SavPhase),intent(inout):: vSavFas(:)
  !
  integer :: I,P,Q,nP
  
  do I=1,size(vSavFas)
    !nP= vMixModel(vSavFas(I)%iModel)%nPole
    nP= vMixModel(I)%nPole
    do P=1,nP
      !--start from compos'n close to end-member P
      vSavFas(i)%tXPole(P,P)= One - 1.0D-3
      do Q=1,nP
        if(Q/=P) vSavFas(i)%tXPole(P,Q)= 1.0D-3/real(nP-1)
      enddo
    enddo
  enddo
  
  return
end subroutine MixPhase_XPole_Init

subroutine Mixture_Minimize( &
& TdgK,Pbar, &
& vFasPur,vMixModel,  &
& vMixFas_Xpole_Init, &
& vMixFas0)
!--
!-- for each mixture phase,
!-- compute the composition X that minimizes G
!-- and compute G at X, to see whether it is -0
!-- (the value istself is useful only for output to trace file)
!--
!-- if there is no phase with negative G,
!-- then the current phase assemblage is stable
!--
  use M_Numeric_Const, only: Ln10
  use M_Safe_Functions,only: FSafe_Exp
  !~ use M_GEM_Vars,   only: TdgK,Pbar
  !
  use M_Optimsolver_Theriak
  use M_MixModel_Optim
  !
  real(dp),        intent(in)   :: TdgK,Pbar
  type(T_Phase),   intent(in)   :: vFasPur(:)
  type(T_MixModel),intent(in)   :: vMixModel(:)
  type(T_SavPhase),intent(inout):: vMixFas_Xpole_Init(:)
  type(T_MixPhase),intent(inout):: vMixFas0(:)
  !
  type(T_MixModel):: MM
  real(dp),allocatable:: vX(:),vXmin(:),vMu(:),vMuMin(:)
  real(dp):: G,Gmin
  integer :: N,I,J,K,P !,iFas
  integer :: Multi
  real(dp):: S
  !
  real(dp):: TolX,DeltaInit
  integer :: its,nCallG
  logical :: OkConverge
  !
  TolX= MixMinim_TolX
  DeltaInit= 0.05D0
  !
  do I=1,size(vMixFas0)
    !
    MM= vMixModel(vMixFas0(I)%iModel)
    N= MM%NPole
    Multi= MM%vMulti(1)
    !
    !print '(A,2I3,1X,A,1X,A)', &
    !& "I,N,FAS,MODEL=", &
    !& I,N,trim(vMixFas0(I)%Name),trim(MM%Name)
    !pause
    !
    allocate(vXmin(N))
    allocate(vMuMin(N))
    !
    if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
      !-------------------------------------------------analytic minimum
      Multi=  MM%vMulti(1)
      !
      S= Zero
      do P=1,N
        vXmin(P)= FSafe_Exp(-vFasPur(MM%vIPole(P))%Grt /real(Multi))
        S= S + vXmin(P)
      enddo
      vXmin(1:N)=  vXmin(1:N) /S
      vMuMin(1:N)= vFasPur(MM%vIPole(1:N))%Grt + Multi*log(vXmin(1:N))
      Gmin= SUM(vXmin(1:N) * vMuMin(1:N))
      !------------------------------------------------/analytic minimum
    else
      !---------------------------------------------numerical minimum(s)
      call MixModel_Optim_SetParams(TdgK,Pbar,MM)
      allocate(Mixmodel_Optim_vMu0rt(N))
      allocate(Mixmodel_Optim_vLPole(N))
      Mixmodel_Optim_vMu0rt(1:N)= vFasPur(MM%vIPole(1:N))%Grt
      Mixmodel_Optim_vLPole(1:N)= .true.
      !
      allocate(vX(N))
      allocate(vMu(N))
      Gmin= 1.0D30
      !
      do J=1,N
        !
        !!! !--- start from compos'n close to end-member J
        !!! vX(J)= One - 1.0D-3
        !!! do K=1,N
        !!!   if(K/=J) vX(K)= 1.0D-3/real(N-1)
        !!! enddo
        !
        vX(1:N)= vMixFas_Xpole_Init(I)%tXPole(J,1:N)
        !
        call Optimsolver_Theriak( & !
        !& MixModel_Optim_GetGibbs,      & !
        !& MixModel_Optim_GetPotentials, & !
        & MixModel_Optim_GetMu, & ! interface
        & FF,                   & ! IN
        & DeltaInit,            & ! IN
        & TolX,                 & ! IN
        & vX,                   & ! INOUT
        & vMu,                  & ! OUT
        & G,                    & ! OUT
        & OkConverge,           & ! OUT
        & its,nCallG)             ! OUT
        !
        if(iDebug>2) write(91,'(I3,2X,A)') its, trim(MM%Name)
        !
        if(WarmRestart) vMixFas_Xpole_Init(I)%tXPole(J,1:N)= vX(1:N)
        !
        if(G<Gmin) then
          Gmin= G
          vXmin(:)= vX(:)
          vMuMin(:)= vMu(:)
        end if
        !
      end do
      !pause
      !
      deallocate(vX)
      deallocate(vMu)
      !
      deallocate(Mixmodel_Optim_vMu0rt)
      deallocate(Mixmodel_Optim_vLPole)
      !--------------------------------------------/numerical minimum(s)
      !
    end if
    !
    if(iDebug>2) write(91,'(A)') "================================="
    !
    vMixFas0(i)%vLPole(1:N)= .true.
    vMixFas0(i)%vXPole(1:N)= vXmin(1:N)
    vMixFas0(i)%Grt= Gmin
    vMixFas0(i)%vLnAct(1:N)= vMuMin(1:N) -vFasPur(MM%vIPole(1:N))%Grt
    !
    deallocate(vXmin)
    deallocate(vMuMin)
    !
  enddo
  
  !------------------------------------------------------------log files
  if(F1>0) then
    write(F1,'(/,A,/)') "Mixture_Minimize,EndMember,Xi,Gi"
    do I=1,size(vMixModel)
      !if(vMixFas0(I)%Grt > 1.D-3) cycle
      MM= vMixModel(I)
      write(F1,'(2A)') "MODEL= ",MM%Name
      do K=1,MM%NPole
        write(F1,'(A,1X,G15.6)') &
        & MM%vNamPole(K),        &
        & vMixFas0(I)%vXPole(K)
      enddo
      write(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
      write(F1,*)
    enddo
    write(F1,*)
  end if
  !
  if(F2>0) then
    do I=1,size(vMixModel)
    !if(vMixFas0(I)%Grt < Zero) then
      MM= vMixModel(I)
      !
      write(F2,'(A,A1)',advance="NO") MM%Name,T_
      do K=1,MM%NPole
        write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%vXPole(K),T_
      enddo
      !do K=1,MM%NPole
      !  write(F2,'(G15.6,A1)',advance="NO") vFasPur(MM%vIPole(K))%Grt,T_
      !enddo
      write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%Grt,T_
      !
    !end if
    enddo
    write(F2,*)
  end if
  !-----------------------------------------------------------/log files
  
  return
end subroutine Mixture_Minimize

subroutine Speciation
end subroutine Speciation

subroutine Write_Log_Entete(vFas)
  use M_IoTools,only: GetUnit
  use M_Files,  only: DirOut
  !
  type(T_Phase),intent(in):: vFas(:)
  !
  call GetUnit(F1)
  open(F1,file=trim(DirOut)//"_gibbs.log")
  write(F1,'(A)') "g of all phases, relative to current assemblage"
  !
  call GetUnit(F2)
  open(F2,file=trim(DirOut)//"_mixcomp.log")
  write(F2,'(A)')   "properties of mixtures with min(g)<0 at each iteration"
  write(F2,'(A,/)') "mix.model.name, Xi, Gmin"
  !
end subroutine Write_Log_Entete

end module M_GEM_Theriak
