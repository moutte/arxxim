module M_Simplex_Theriak_Grid
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  use M_Kinds
  use M_Trace,     only: iDebug,fTrc,T_,Stop_,Pause_
  !
  use M_T_Component,  only: T_Component
  use M_T_Species, only: T_Species
  use M_T_Phase,   only: T_Phase
  use M_T_MixPhase,only: T_MixPhase
  use M_T_MixModel,only: T_MixModel,Mix_Molecular
  use M_T_MixModel,only: MaxPole
  use M_GEM_Vars,  only: T_SavPhase
  !
  implicit none
  !
  private
  !
  public:: Simplex_Theriak_Grid
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
  !----------------------------------------------------private variables
  integer:: F1,F2,FF
  !
  type(T_SavPhase):: SavPhaseZero
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
  logical,parameter:: WarmRestart= .false.  !!.true.  !!
  !--------------------------------------------------------------------/
  !
contains

subroutine Simplex_Theriak_Grid
!--
!-- construction of a 2D phase diagram
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  use M_IoTools,     only: GetUnit
  use M_Dtb_Const,   only: T_CK
  ! use M_Simplex_Vars,only: Simplex_Vars_Alloc,Simplex_Vars_Clean
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
  ! use M_Simplex_Vars,only: tSimplex
  !
  !--for Global_TP_Update
  use M_Global_Tools,only: Global_TP_Update
  use M_Global_Vars, only: vSpcDtb,vSpc,vMixModel
  use M_Global_Vars, only: vDiscretModel,vDiscretParam,vMixFas
  !---------------------------------------------------------------------
  !
  type:: T_Node2D
    real(dp):: X,Y
    integer :: Nord,Sud,Est,West
  end type T_Node2D
  !
  type:: T_Lis_Node2D
    type(T_Node2D)::Value
    type(T_Lis_Node2D),pointer::Next
  end type T_Lis_Node2D
  type(T_Lis_Node2D),pointer:: Lis_Node2D
  !
  type(T_Node2D),allocatable:: vNode2D(:)
  integer :: nNodX, nNodY
  !
  ! type(T_Lis_XY_Node),pointer:: LisCur, LisPrev
  !
  ! if(N==1) then
  !   allocate(Lis_XY_Node)
  !   nullify(Lis_XY_Node%next)
  !   Lis_XY_Node%Value= XY_Node
  !   LisCur=> Lis_XY_Node
  ! else
  !   allocate(LisCur%next)
  !   nullify(LisCur%next%next)
  !   LisCur%next%Value= XY_Node
  !   LisCur=> LisCur%next
  ! end if
  !
  integer :: I,J,K,iPath
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
  
  if(idebug>1) write(fTrc,'(/,A)') "< Simplex_Theriak_Grid"
  !
  nNodX= 20
  nNodY= 20
  allocate(vNode2D(nNodX*nNodY))
  TdgK= TdgK0
  Pbar= Pbar0
  K= 0
  do I=1,nNodX
  end do
  !
  nC= size(vCpnGEM)
  nFpur= size(vFas) ! + size(vMixModel)
  nMix= size(vMixModel)
  !
  F1= 0
  F2= 0
  FF= 0
  !
  !---------------------------------read path parameters from PATH block
  ! ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? ??? 
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
    SavPhaseZero%iModel=      0
    SavPhaseZero%nFas=        0
    SavPhaseZero%tXPole(:,:)= Zero
    SavPhaseZero%vMole(:)=    Zero
    SavPhaseZero%vGrt0(:)=    Zero
    SavPhaseZero%vVol0(:)=    Zero
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

  do
    !
    TdgK0= TdgK ! previous T,P values
    Pbar0= Pbar ! 
    !
    !----------------------if change T,P update T,P-dependent properties
    ! TdgK= vTPpath(iPath)%TdgC +T_CK
    ! Pbar= vTPpath(iPath)%Pbar
    !
    if(TdgK0 /= TdgK .or. Pbar0 /= Pbar) &
    & call Global_TP_Update( &
    & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
    & vSpc,vMixModel,vMixFas,vFas)
    !------------------------------------------------------------------/
    !
    !-------------------------------------------------system composition
    !change ... ! do J=1,nC
    !change ... !   if(vLPath(J)) then
    !change ... !     vCpnGEM(J)%Mole= tPathData(J,iPath)
    !change ... !   end if
    !change ... ! end do
    !------------------------------------------------------------------/
    !
    !change ... ! tResult(1,     iPath)= TdgK -T_CK
    !change ... ! tResult(2,     iPath)= Pbar
    !change ... ! tResult(3:nC+2,iPath)= vCpnGEM(1:nC)%Mole
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
          !print *,vSavPhase(i)%vVol0(1)
          if(vSavPhase(i)%nFas > 1) call Check_vXMean(i)
        end do
        !pause !!
      end if
      !
      !call Path_StoreResult(iPath,vCpnGEM(:)%Mole,TdgK,Pbar)
      !
    end if
    !-----------------------/store results in tables tResult, tResultMix
    !
  end do
  !
  call WritePhases(       &
  & vSpc,vFas,vMixModel,  &
  & vCpnGEM,              &
  & dimPath,vSimplex_Ok)
  !
  if(size(vMixModel)>0) &
  & call WriteMixtures( &
  & vSpc,vFas,vMixModel,    &
  & vCpnGEM,                &
  & DimPath,                &
  & vSimplex_Ok)
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
  if(idebug>1) write(fTrc,'(A,/)') "</ Simplex_Theriak_Path"
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
  end do
  vX(:)= vX(:) /sum(vX(:))
  !
  G= MixPhase_Grt( &
  & TdgK,Pbar,    &
  & nP, &
  & vMixModel(J), &
  & vSavPhase(i)%vGrt0(1:nP), &
  & vX(1:nP))
  !
  deallocate(vX)

  !! print *,"G(vXMean)..",vMixModel(J)%Name, G  ;  pause

  return
end subroutine Check_vXMean

end subroutine Simplex_Theriak_Grid

subroutine Simplex_GEM(iError)
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  use M_Simplex_Calc,only: Simplex_Calc
  !! use M_Simplex_Vars,only: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  use M_Global_Vars, only: vFas,vMixModel
  use M_GEM_Vars,    only: vCpnGEM,tStoikioGEM
  use M_GEM_Vars,    only: TdgK,Pbar
  !
  !! use M_Simplex_Vars,only: iPosV,tSimplex
  !
  integer, intent(out):: iError
  !
  real(dp),        allocatable:: tStoikio(:,:)
  type(T_Phase),   allocatable:: vFasPur(:)
  type(T_GemPhase),allocatable:: vFas0(:)
  type(T_MixPhase),allocatable:: vMixFas0(:)
  type(T_SavPhase),allocatable:: vMixFas_Xpole_Init(:)
  !
  ! type(T_MixModel):: MM
  !
  real(dp),allocatable:: tSimplex(:,:)
  integer, allocatable:: IZROV(:),IPOSV(:)
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
  nC=    size(vCpnGEM)
  nFpur= size(vFas) !-- at this point, only pure phases are in vFas ==
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
  end do
  !
  if(nMix>0) then
    do K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    end do
    call MixPhase_XPole_Init(vMixModel,vMixFas_Xpole_Init)
  end if

  !------------------------------------------------------initial simplex
  nFmixSpl= 0
  nF= nFpur
  allocate(tSimplex(0:nC+1,0:nF))  ;  tSimplex=Zero
  allocate(IPOSV(1:nC))            ;  IPOSV= 0
  allocate(IZROV(1:nF))            ;  IZROV= 0
  !
  !! call Simplex_Vars_Alloc(nC,nF)
  !
  tSimplex(1:nC, 0)=  vCpnGEM(1:nC)%Mole
  tSimplex(0,    1:nF)= -vFas0(1:nF)%Grt
  tSimplex(1:nC, 1:nF)= -transpose(tStoikio(1:nF,1:nC))
  !
  !! call Simplex_Calc(iError)
  call Simplex_Calc( &
  & nC,nF,           &
  & tSimplex,        &
  & IZROV,IPOSV,     &
  & iError) !!,n1,n2)
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
      if(all(vMixFas0(:)%Grt >=-G_Iota)) then
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
      end do
      !
      ! call Simplex_Vars_Clean
      if(allocated(tSimplex))  deallocate(tSimplex)
      if(allocated(IZROV))     deallocate(IZROV)
      if(allocated(IPOSV))     deallocate(IPOSV)
      ! call Simplex_Vars_Alloc(nC,nF)
      allocate(tSimplex(0:nC+1,0:nF))  ;  tSimplex=Zero
      allocate(IPOSV(1:nC))            ;  IPOSV= 0
      allocate(IZROV(1:nF))            ;  IZROV= 0
      !
      tSimplex(1:nC,0   )=  vCpnGEM(1:nC)%Mole
      tSimplex(0,   1:nF)= -vFas0(1:nF)%Grt
      tSimplex(1:nC,1:nF)= -transpose(tStoikio(1:nF,1:nC))
      !
      !! call Simplex_Calc(iError)
      call Simplex_Calc( &
      & nC,nF,           &
      & tSimplex,        &
      & IZROV,IPOSV,     &
      & iError) !!,n1,n2)
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
        end do
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
    end do
    !
    if(iDebug>2) print '(A,I3)',"Iter= ",Iter
    !
  end if
  !
  call GEM_SaveResults(nFpur,Iter,vMixModel)
  !
  !--------------------------------------------------------DesAllocation
  ! call Simplex_Vars_Clean
  if(allocated(tSimplex))  deallocate(tSimplex)
  if(allocated(IZROV))     deallocate(IZROV)
  if(allocated(IPOSV))     deallocate(IPOSV)
  !
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
    end do
    
    write(6,'(A)') "MIX="
    do K=1,nC
      I= IPOSV(K)
      if(vFas0(I)%iModel/=0) then
        write(6,'(4X,2I3,2X,G15.6,A)',advance="NO") &
        & I,K,tSimplex(K,0)," X="
        J= vFas0(I)%iModel
        do P=1,vMixModel(J)%NPole
          write(6,'(F7.3,1X)',advance="NO") vFas0(I)%vXPole(P)
        end do
        write(6,'(2X,A)') trim(vFas0(I)%NamFs)
      end if
    end do
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
    end do

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
        end do
        !
      end if
      !
    end do

    !--------------------------------------------------------------trace
    if(iDebug>2) then
      do K=1,nC
        I= IPOSV(K)
        if(vFas0(I)%iModel/=0) then
          write(21,'(A,A1,G15.6,A1)',advance="NO") &
          & trim(vFas0(I)%NamFs),T_, &
          & tSimplex(K,0),       T_
          J= vFas0(I)%iModel
          do P=1,vMixModel(J)%NPole
            write(21,'(G15.6,A1)',advance="NO") vFas0(I)%vXPole(P),T_
          end do
        end if
      end do
      write(21,*)
    end if
    !-------------------------------------------------------------/trace

    return
  end subroutine StoreResults

end subroutine Simplex_GEM

real(dp) function MixPhase_Grt(TdgK,Pbar,nP,MM,vMu0rt,vX)
  use M_T_MixModel,only: MixModel_Activities

  real(dp),        intent(in):: TdgK,Pbar
  integer,         intent(in):: nP
  type(T_MixModel),intent(in):: MM        ! mixing model
  real(dp),        intent(in):: vMu0rt(:) !
  real(dp),        intent(in):: vX(:)     ! phase composition
  !
  real(dp):: vLGam(nP),vLIdeal(nP),vLnAct(nP)
  logical :: vLPole(nP)
  real(dp):: G
  integer :: i
  logical :: Ok
  character(len=30):: Msg

  vLPole(:)= (vX(:)>Zero)

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
  end do
  
  !--- mixtures --
  vFasMolMix(:)= Zero
  do I=1,size(vSavPhase)
    if(vSavPhase(I)%nFas >0) then
      vFasIsPresent(nFpur +I)= .true.
      do J=1,vSavPhase(I)%nFas
        vFasMolMix(I)= vFasMolMix(I) +vSavPhase(I)%vMole(J)
      end do
    end if
  end do
  
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
  end do
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
        & sum( vFas(MM%vIPole(1:N))%vStoik(iCp) &
        &    * vMixFas(iMix)%vXPole(1:N) )
      end do
      !
      if(iDebug>2) print *,"ADDED= ",trim(vFas(iMF)%NamFs)
      !
    end if
  
  end do
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
  !use M_Simplex_Vars, only: IPOSV
  !
  integer, intent(in):: nC,nF
  integer, intent(in):: IPOSV(:)
  type(T_GemPhase),intent(inout):: vFas0(:)
  !
  real(dp),allocatable:: tTransform(:,:)
  integer, allocatable:: vIndx(:)
  real(dp),allocatable:: vY(:)
  real(dp),allocatable:: vGrt(:),vGrt0(:)
  !
  integer :: I,K
  logical :: bSingul
  real(dp):: D
  !
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
  end do
  !--------------------------------------------------------------------/
  !
  call LU_Decomp(tTransform,vIndx,D,bSingul)
  if(bSingul) call Stop_("SINGUL IN Gibbs_Change")
  !
  do I=1,nF
    vY(1:nC)= vFas0(I)%vStoik(1:nC)
    !-----------------------Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
    call LU_BakSub(tTransform,vIndx,vY)
    vGrt(I)= vFas0(I)%Grt - sum(vY(:)*vGrt0(:))
  end do
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
      end do
    end do
  end do
  
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
  !! use M_GEM_Vars,   only: TdgK,Pbar
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
      end do
      vXmin(1:N)=  vXmin(1:N) /S
      vMuMin(1:N)= vFasPur(MM%vIPole(1:N))%Grt + Multi*log(vXmin(1:N))
      Gmin= sum(vXmin(1:N) * vMuMin(1:N))
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
        !! !--- start from compos'n close to end-member J
        !! vX(J)= One - 1.0D-3
        !! do K=1,N
          !! if(K/=J) vX(K)= 1.0D-3/real(N-1)
        !! end do
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
        write(77,'(I3,2X,A)') its, trim(MM%Name)
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
    write(77,'(A)') "================================="
    !
    vMixFas0(i)%vLPole(1:N)= .true.
    vMixFas0(i)%vXPole(1:N)= vXmin(1:N)
    vMixFas0(i)%Grt= Gmin
    vMixFas0(i)%vLnAct(1:N)= vMuMin(1:N) -vFasPur(MM%vIPole(1:N))%Grt
    !
    deallocate(vXmin)
    deallocate(vMuMin)
    !
  end do
  
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
      end do
      write(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
      write(F1,*)
    end do
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
      end do
      !do K=1,MM%NPole
      !  write(F2,'(G15.6,A1)',advance="NO") vFasPur(MM%vIPole(K))%Grt,T_
      !end do
      write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%Grt,T_
      !
    !end if
    end do
    write(F2,*)
  end if
  !-----------------------------------------------------------/log files
  
  return
end subroutine Mixture_Minimize

subroutine WritePhases( &
& vSpc,vFas,vMixModel,  &
& vCpnGEM,              &
& DimPath,              &
& vSimplex_Ok)
!--
!-- write tabulated results for all phases found stable at any step
!--
  use M_Files,        only: DirOut
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: GetUnit
  !use M_Numeric_Tools,only: iMaxLoc_R
  !
  !use M_Global_Vars,only: vSpc,vFas,vMixModel
  !use M_GEM_Vars,   only: vCpnGEM
  !---------------------------------------------------------------------
  type(T_Species),  intent(in):: vSpc(:)
  type(T_Phase),    intent(in):: vFas(:)
  type(T_MixmOdel), intent(in):: vMixModel(:)
  type(T_Component),intent(in):: vCpnGEM(:)
  !
  integer, intent(in):: DimPath
  logical, intent(in):: vSimplex_Ok(DimPath)
  !---------------------------------------------------------------------
  integer :: iPath,iFs,I,J,K
  integer :: nC,nF,nP,nFpur,nFmix
  integer :: F
  !real(dp):: vX(MaxPole)
  real(dp):: Tot,X,Y
  type(T_SavPhase):: Fas0
  type(T_MixModel):: MM
  !
  logical,parameter:: ComputeXMean= .true. !.false.
  !---------------------------------------------------------------------
  nC= size(vCpnGEM)
  nF= size(vFasIsPresent)
  nFpur= size(vFas)
  nFmix= size(vMixModel)
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_phase_mole.restab")
  call  Write_Title(F)
  !
  call GetUnit(F1)
  open(F1,file=trim(DirOut)//"_phase_cm3.restab")
  call  Write_Title(F1)
  !
  call GetUnit(F2)
  open(F2,file=trim(DirOut)//"_phase_grams.restab")
  call  Write_Title(F2)
  !
  K=0
  do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then
      !
      K=K+1
      !
      write(F, '(I3,A1)',      advance="no") K,T_  != count
      write(F1,'(I3,A1)',      advance="no") K,T_  != count
      write(F2,'(I3,A1)',      advance="no") K,T_  != count
      !
      write(F,'(2(F12.3,A1))',advance="no") &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      do iFs=1,nC != (3:nC+2)                     != components
        write(F,'(G15.6,A1)',advance="no") tResult(2+iFs,iPath),T_
      end do
      !
      write(F1,'(2(F12.3,A1))',advance="no") &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      do iFs=1,nC != (3:nC+2)                     != components
        write(F1,'(G15.6,A1)',advance="no") tResult(2+iFs,iPath),T_
      end do
      !
      write(F2,'(2(F12.3,A1))',advance="no") &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      do iFs=1,nC != (3:nC+2)                     != components
        write(F2,'(G15.6,A1)',advance="no") tResult(2+iFs,iPath),T_
      end do
      !
      !! Tot=sum(tResult(1:nC,iPath))
      !! do iFs=1,nC
        !! write(F,'(G15.6,A1)',advance="no") tResult(iFs,iPath)/Tot,T_
      !! end do
      
      !-----------------------------------------------------mole numbers
      do iFs=1,nF != (nC+3:nC+nF+2)
        if(vFasIsPresent(iFs)) &
        & write(F,'(G15.3,A1)',advance="no") &
        ! nr'moles phase          *nr'oxygen in formula
        ! & tResult(nC+iFs,iPath)*tFormula(iFs,1),T_
        & tResult(2+nC+iFs,iPath),T_
      end do
      write(F,*)
      !----------------------------------------------------/mole numbers

      !---------------------------------------------volume/cm3, weight/g
      do iFs=1,nFpur != (nC+3:nC+nF+2)
        if(vFasIsPresent(iFs)) &
        & write(F1,'(G15.3,A1)',advance="no") &
        & tResult(2+nC+iFs,iPath)*vFas(iFs)%VolM3*1.D6,T_
        if(vFasIsPresent(iFs)) &
        & write(F2,'(G15.3,A1)',advance="no") &
        & tResult(2+nC+iFs,iPath)*vFas(iFs)%WeitKg*1.D3,T_
      end do
      do iFs=1,nFmix
        if(vFasIsPresent(nFpur +iFs)) then
        
          Fas0= tResultMix(iFs,iPath)
          MM= vMixModel(Fas0%iModel)
          X= 0.D0
          Y= 0.D0
          do J=1,MM%nPole
            do I=1,MM%nPole
              X= X + vSpc(MM%vIPole(I))%V0     * Fas0%tXPole(J,I)
              Y= Y + vSpc(MM%vIPole(I))%WeitKg * Fas0%tXPole(J,I)
            end do
          end do
          
          !X= dot_product( &
          !& tResultMix(iFs,iPath)%tXPole(1,:), &
          !& tResultMix(iFs,iPath)%vVol0(:))
          
          write(F1,'(G15.3,A1)',advance="no") &
          & tResult(2+nC+nFpur+iFs,iPath)*X*1.D6,T_
          write(F2,'(G15.3,A1)',advance="no") &
          & tResult(2+nC+nFpur+iFs,iPath)*Y*1.D3,T_
          
        end if
      end do
      write(F1,*)
      write(F2,*)
      !--------------------------------------------/volume cm3, weight/g
    end if !!if(vSimplex_Ok(iPath))

  end do !!do iPath
  !
  write(F,'(A1)') "_"
  !
  close(F)
  close(F1)
  close(F2)
  !
  if(iDebug>2) &
  & print '(2A)', "Phases: Results in ",trim(DirOut)//"_phase_xx.restab"
  !
contains

  subroutine Write_Title(FF)
    integer,intent(in):: FF
    !--------------------------------------------------------title lines
    write(FF,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
    !
    !--component names
    do i=1,nC
      write(FF,'(A,A1)',advance="no") trim(vCpnGEM(i)%NamCp)//"_cpn", T_
    end do
    !
    !--pure phase names
    do iFs=1,nFpur
      if(vFasIsPresent(iFs)) &
      & write(FF,'(A,A1)',advance="no") trim(vFas(iFs)%NamFs), T_
    end do
    !--mixture phase names
    do iFs=1,nFmix
      if(vFasIsPresent(nFpur +iFs)) &
      & write(FF,'(A,A1)',advance="no") trim(vMixModel(iFs)%Name), T_
    end do
    write(FF,*)
    !-------------------------------------------------------/title lines
  end subroutine Write_Title

end subroutine WritePhases

subroutine WriteMixtures( &
& vSpc,vFas,vMixModel,    &
& vCpnGEM,                &
& DimPath,                &
& vSimplex_Ok)
!--
!-- write tabulated results for all phases found stable at any step
!--
  use M_Files,        only: DirOut
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: GetUnit
  use M_Numeric_Tools,only: iMaxLoc_R
  !
  !---------------------------------------------------------------------
  type(T_Species),  intent(in):: vSpc(:)
  type(T_Phase),    intent(in):: vFas(:)
  type(T_MixmOdel), intent(in):: vMixModel(:)
  type(T_Component),intent(in):: vCpnGEM(:)
  integer,          intent(in):: DimPath
  logical,          intent(in):: vSimplex_Ok(DimPath)
  !---------------------------------------------------------------------
  integer :: iPath,iFs,I,J,K,P,Q
  integer :: nP,nFmix,nFPur,nPP
  integer :: FMIX
  real(dp):: vX(MaxPole)
  real(dp):: Tot
  type(T_SavPhase):: Fas0,Fas1
  type(T_MixModel):: MM
  !
  logical,parameter:: ComputeXMean= .true. !.false.
  
  call GetUnit(FMIX)
  open(FMIX,file=trim(DirOut)//"_mixtures.restab")
  !
  nFpur= size(vFas)
  nFmix= size(vMixModel)
  !
  !-----------------------------------------------------------title line
  write(FMIX,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
  !
  do iFs=1,nFmix
    !
    if(vFasIsPresent(nFpur +iFs)) then
      !
      !Fas= tResultMix(iFs,iPath)
      MM= vMixModel(iFs)
      nP= MM%nPole
      !if(Fas0%iModel /= 0) then
      if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
        nPP= 1
      else
        nPP= MM%nPole
      end if
      do I=1,nPP
        write(FMIX,'(A,A1)',advance="NO") trim(MM%Name),T_
        do P=1,nP
          write(FMIX,'(A,A1)',advance="NO") &
          & trim(vSpc(MM%vIPole(P))%NamSp),T_
        end do
      end do
      !
    end if
  
  end do
  write(FMIX,*)
  !----------------------------------------------------------/title line
  
  K=0
  do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then
      !
      K=K+1
      !------------------------------compute mean comp'n of each mixture
      !----------------------------------------------and its free energy
      if(ComputeXMean) then
        !
        do J=1,nFmix
          !
          Fas0=        tResultMix(J,iPath)
          Fas1=        SavPhaseZero
          Fas1%iModel= Fas0%iModel
          nP=          vMixModel(Fas0%iModel)%nPole
          Fas1%nFas=   nP
          !
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          !if(Fas0%iModel /= 0) then
          do P=1,Fas0%nFas
            do Q=1,Fas1%nFas

              vX(1:nP)= abs( Fas0%tXPole(P,1:nP) - Fas0%tXPole(Q,1:nP) )

              if(maxval(vX(1:nP)) < MixMinim_TolX *1.0D2) then
              
                Tot= Fas0%vMole(P) +Fas1%vMole(Q)
                vX(1:nP)= Fas0%vMole(P) *Fas0%tXPole(P,1:nP) &
                &       + Fas1%vMole(Q) *Fas0%tXPole(Q,1:nP)
                Fas1%tXPole(Q,1:nP)= vX(1:nP) /Tot
                Fas1%vMole(Q)= Tot
                
                exit
              
              end if

            end do
          end do
          !
          tResultMix(J,iPath)= Fas1
          !end if
          !
        end do
        !
      end if
      !----------------------------------/compute mean comp'n of mixture

      !---------------------------------------write mixture compositions
      !-----------------------------------------------------in file FMIX
      write(FMIX,'(I3,A1)',     advance="no") K,T_   != count
      !
      write(FMIX,'(2(F7.2,A1))',advance="no") &      !
      & tResult(1,iPath), T_, &                      != TdgC
      & tResult(2,iPath), T_                         != Pbar

      do J=1,nFmix
        if(vFasIsPresent(nFpur +J)) then
          Fas0= tResultMix(J,iPath)
          MM= vMixModel(Fas0%iModel)
          nP= MM%nPole
          !if(Fas0%iModel /= 0) then
          if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
            nPP= 1
          else
            nPP= MM%nPole
          end if
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          do I=1,nPP
            write(FMIX,'(G15.6,A1)',advance="NO") Fas0%vMole(I),T_
            do P=1,nP
              write(FMIX,'(F7.3,A1)',advance="NO") &
              & Fas0%tXPole(I,P),T_
            end do
          end do
          !end if
        end if
      end do
      !
      write(FMIX,*)
      !--------------------------------------/write mixture compositions
      !
    end if !!if(vSimplex_Ok(iPath))

  end do !!do iPath
  !
  if(iDebug>2) &
  & print '(2A)', "Mixtures: Results in ",trim(DirOut)//"_mixtures.restab"
  !
end subroutine WriteMixtures

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

end module M_Simplex_Theriak_Grid
