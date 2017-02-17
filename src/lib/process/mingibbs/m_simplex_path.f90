module M_Simplex_Path
!-- series of equilibrium calculations on pure phases using simplex method
!-- under variable T,P conditions, or under variable bulk composition
  !
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  !
  implicit none
  !
  private
  !
  public:: Simplex_Path
  ! series of equilibrium calculations ON pure PHASES using simplex method
  !
  real(dp),allocatable,public:: tSimplexResult(:,:)
  logical, allocatable,public:: vSimplex_Ok(:)
  logical, allocatable,public:: vFasIsPresent(:)
  !
contains

subroutine Simplex_Path(Cod)
!--
!-- series of equilibrium calculations on pure phases
!-- using simplex method
!-- under variable T,P conditions,
!-- or under variable bulk composition
!--
  use M_IoTools,     only: OutStrVec
  use M_Files,       only: NamFInn
  use M_Dtb_Const,   only: T_CK, R_JK
  use M_Global_Tools,only: Global_TP_Update
  use M_TPcond_Read, only: TPpath_Read
  use M_Path_Read,   only: Path_ReadMode, Path_ReadParam_new
  use M_Path_Vars,   only: Path_Vars_Clean
  !
  use M_Simplex_Build,only: Simplex_CloseFiles
  use M_Simplex_Calc, only: Simplex_Calc
  !
  use M_Global_Vars, only: vFas,vSpcDtb,vSpc,vMixModel
  use M_Global_Vars, only: vDiscretModel,vDiscretParam,vMixFas
  !
  use M_Path_Vars,   only: vTPpath,vLPath,tPathData,DimPath
  !
  use M_GEM_Vars,    only: vCpnGEM,TdgK,Pbar,tStoikioGEM
  !---------------------------------------------------------------------
  character(len=*),intent(in) :: Cod
  !---------------------------------------------------------------------
  integer :: I,J,nC,nF
  real(dp):: TdgK0,Pbar0
  logical :: Ok
  character(len=3) :: PathMod3
  character(len=80):: Msg
  !
  real(dp),allocatable:: tSimplex(:,:)
  integer, allocatable:: IZROV(:),IPOSV(:)
  !---------------------------------------------------------------------
  Ok= .true.
  
  nC= size(vCpnGEM)
  nF= size(vFas)
  
  !! call Simplex_Vars_Alloc(nC,nF)
  if(allocated(tSimplex))  deallocate(tSimplex)
  if(allocated(IZROV))     deallocate(IZROV)
  if(allocated(IPOSV))     deallocate(IPOSV)
  !
  allocate(tSimplex(0:nC+1,0:nF))  ;  tSimplex=Zero
  allocate(IPOSV(1:nC))            ;  IPOSV= 0
  allocate(IZROV(1:nF))            ;  IZROV= 0
  
  ! if(iDebug>2) then
  ! print *,'done Simplex_Vars_Alloc' !  ;  call pause_
  ! end if
  
  ! tSimplex=
  ! row 0 =    Gibbs energy of phases 1:nF at T,P
  ! column 0 = bulk compos'n
  ! main =     stoikiometry matrix
  tSimplex(1:nC,1:nF)= -transpose(tStoikioGEM(1:nF,1:nC)) !stoikio
  tSimplex(0,   1:nF)= -vFas(1:nF)%Grt  !Gibbs energy
  tSimplex(1:nC,0   )=  vCpnGEM(1:nC)%Mole !bulk compos'n
  
  if(allocated(vFasIsPresent)) deallocate(vFasIsPresent)
  allocate(vFasIsPresent(1:nF))  ; vFasIsPresent=.false.
  
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  
  select case(trim(Cod))
  
  case default
    Msg= trim(Cod)//"= invalid code in SIMPLEX PATH"
    Ok= .false.
    print *,trim(Msg)
    return
  
  case("PATH")
  !------------------------------------------------------------case PATH
  !--- calculation at fixed T,P,
  !--- changing step by step the amount of components' mole numbers
    !
    !--------------------------- read path parameters from PATH block --
    call Path_ReadMode(NamFInn,PathMod3,Ok,Msg)
    !
    if(PathMod3 /= "CHG") then
      Msg= "Global equilibrium paths: only in CHANGE mode"
      Ok= .false.
      print *,trim(Msg)
      return
    end if
    !
    call Path_ReadParam_new( &
    & NamFInn,   &
    & PathMod3,  &
    & vCpnGEM,   &
    & TdgK,Pbar, &
    & Ok,Msg)
    !
    if(iDebug>2) then
    print *,'done Path_ReadParam'  ;  call pause_
    end if
    !
    allocate(vSimplex_Ok(1:dimPath))  ; vSimplex_Ok=.false.
    !
    if(allocated(tSimplexResult)) deallocate(tSimplexResult)
    allocate(tSimplexResult(1:nC+nF+2,1:dimPath)); tSimplexResult=Zero
    !
    !~ call Global_TP_Update( &
    !~ & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
    !~ & vSpc,vMixModel,vMixFas,vFas)
    !
    do i=1,dimPath
      !
      TdgK0= TdgK
      Pbar0= Pbar
      !--------------------if change T,P update T,P-dependent properties
      TdgK= vTPpath(I)%TdgC +T_CK
      Pbar= vTPpath(I)%Pbar
      !
      if(TdgK0 /= TdgK .or. Pbar0 /= Pbar) &
      & call Global_TP_Update( &
      & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
      & vSpc,vMixModel,vMixFas,vFas)
      !
      do J=1,size(vSpc)
        write(31,'(A,A1,2(F7.1,A1),2(G15.7,A1))') &
        & trim(vSpc(J)%NamSp),T_, &
        & TdgK-T_CK,          T_, &
        & Pbar,               T_, &
        & vSpc(J)%G0rt,       T_, &
        & vSpc(J)%v0,         T_ 
      end do
      !----------------------------------------------------------------/
      !
      !--- system composition --
      do J=1,nC
        if(vLPath(J)) then
          vCpnGEM(J)%Mole= tPathData(J,I)
        end if
      end do
      !---/system composition --
      !
      !--------------------------------------------rebuild tSimplex(:,:)
      ! main = stoikiometry matrix
      tSimplex(1:nC,1:nF)= -transpose(tStoikioGEM(1:nF,1:nC))
      ! first row= Gibbs energy of phases 1:nF at 'point' iTP
      tSimplex(0,1:nF)= -vFas(1:nF)%Grt
      ! first column= bulk compos'n
      tSimplex(1:nC,0)=  vCpnGEM(1:nC)%Mole
      !----------------------------------------------------------------/
      
      if(iDebug==1) print *,"Simplex_Run_Step=",i
      if(iDebug>2) then
      print *,'call Simplex_Run'  ;  call pause_
      end if
      
      !--------------------------------------------------------------run
      call Simplex_Run( &
      & i,              &
      & nC,nF,          &
      & TdgK,Pbar,      &
      & vFas,           &
      & vCpnGEM,        &
      & tStoikioGEM,    &
      & tSimplex,       &
      & IZROV,IPOSV,    &
      & vSimplex_Ok(i))
      !----------------------------------------------------------------/
      !
    end do
  !-----------------------------------------------------------/case PATH
  !
  case("TP")
  !--------------------------------------------------------------case TP
  !--- calculation at fixed composition,
  !--- changing the T,P condition along the TP.TABLE
    !
    !! call Dtb_Tabulate("GIBBS")
    !
    if(iDebug>1) call Trace_1
    !
    call TPpath_Read
    !
    dimPath= size(vTPpath)
    !
    if(dimPath<2) then
      print '(A)',"NO TP.TABLE available for Simplex TP.run"
      return
    end if
    !
    allocate(vSimplex_Ok(1:dimPath))
    vSimplex_Ok=.false.
    !
    if(allocated(tSimplexResult)) deallocate(tSimplexResult)
    allocate(tSimplexResult(1:nC+nF+2,1:dimPath))
    tSimplexResult=Zero
    !
    do i=1,dimPath
      !
      !---------------------change T,P ; update T,P-dependent properties
      TdgK= vTPpath(I)%TdgC +T_CK
      Pbar= vTPpath(I)%Pbar
      !
      call Global_TP_Update( &
      & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
      & vSpc,vMixModel,vMixFas,vFas)
      !----------------------------------------------------------------/
      !
      !--------------------------------------------rebuild tSimplex(:,:)
      !-------------------------------------- main - stoikiometry matrix
      tSimplex(1:nC,1:nF)= -transpose(tStoikioGEM(1:nF,1:nC))
      !------------------------------------- first column- bulk compos'n
      tSimplex(1:nC,0)=  vCpnGEM(1:nC)%Mole
      !----------- first row- Gibbs energy of phases 1:nF at 'point' iTP
      tSimplex(0,1:nF)= -vFas(1:nF)%Grt
      !
      if(iDebug==1) print *,"Simplex_Run_Step=",i
      if(iDebug==4) then
        write(fTrc,'(I3,A1,2(G15.6,A1))',advance="NO") &
        & I,T_, TdgK-T_CK,T_, Pbar,T_
        call OutStrVec(fTrc,vFas%Grt)
      end if
      !-------------------------------------------/rebuild tSimplex(:,:)
      !
      !--------------------------------------------------------------run
      call Simplex_Run( &
      & i,              &
      & nC,nF,          &
      & TdgK,Pbar,      &
      & vFas,           &
      & vCpnGEM,        &
      & tStoikioGEM,    &
      & tSimplex,       &
      & IZROV,IPOSV,    &
      & vSimplex_Ok(i))
      !old! call Simplex_Run( &
      !old! & i,TdgK,Pbar, &
      !old! & vFasIsPresent,tSimplexResult, &
      !old! & vSimplex_Ok(i))
      !-------------------------------------------------------------/run
      !
    end do
    !
    deallocate(vTPpath)
  !-------------------------------------------------------------/case TP
  end select
  !
  !call Files_Close
  !
  call WriteSysComp(dimPath,vCpnGEM) !,vSimplex_Ok,vFasIsPresent,tSimplexResult)
  !
  if(allocated(vSimplex_Ok)) deallocate(vSimplex_Ok)
  deallocate(tSimplexResult)
  deallocate(vFasIsPresent)
  !
  if(allocated(tStoikioGEM)) deallocate(tStoikioGEM)
  !
  if(allocated(vTPpath)) deallocate(vTPpath)
  !
  !call Simplex_Vars_Clean
  if(allocated(tSimplex))  deallocate(tSimplex)
  if(allocated(IZROV))     deallocate(IZROV)
  if(allocated(IPOSV))     deallocate(IPOSV)
  !
  call Simplex_CloseFiles
  !
  call Path_Vars_Clean
  !
contains

subroutine Trace_1

  write(fTrc,'(/,A)') "SPL2"
  write(fTrc,'(3(A,A1))',advance="NO") "Count",T_,"TdgC",T_,"Pbar",T_
  do I=1,size(vFas)
    write(fTrc,'(A15,A1)',advance="NO") vFas(I)%NamFs,T_
  end do
  write(fTrc,*)

  do I=1,size(vFas) !! check phase / species matching ...!!
    write(fTrc,'(A15,A1)',advance="NO") vSpc(vFas(I)%iSpc)%NamSp,T_
  end do
  write(fTrc,*)

end subroutine Trace_1

end subroutine Simplex_Path

subroutine Simplex_Run( &
& iCount,               &
& nC,nF,                &
& TdgK,Pbar,            &
& vFas,                 &
& vCpn,                 &
& tStoikioCpn,          &
& tSimplex,             &
& IZROV,IPOSV,          &
& Ok)

  use M_T_Phase,     only: T_Phase
  use M_T_Component, only: T_Component
  use M_Simplex_Calc,only: Simplex_Calc
  !---------------------------------------------------------------------
  integer,          intent(in) :: iCount
  integer,          intent(in) :: nC,nF
  real(dp),         intent(in) :: TdgK,Pbar
  type(T_Phase),    intent(in) :: vFas(:)
  type(T_Component),intent(in) :: vCpn(nC)
  real(dp),         intent(in) :: tStoikioCpn(nF,nC)
  real(dp),         intent(inout) :: tSimplex(0:nC+1,0:nF)
  integer,          intent(inout) :: IZROV(1:nF)
  integer,          intent(inout) :: IPOSV(1:nC)
  logical,          intent(out):: Ok
  !---------------------------------------------------------------------
  integer :: iError
  !---------------------------------------------------------------------
  
  !! call Simplex_Calc(iError) ! &
  call Simplex_Calc( &
  & nC,nF,           &
  & tSimplex,        &
  & IZROV,IPOSV,     &
  & iError) !!,n1,n2)
  !
  Ok= iError==0

  select case(iError)

  case(0)

    ! call Simplex_WriteResult(TdgK,Pbar,vCpn,tStoikio)
    
    call Simplex_WriteResult( &
    & TdgK,Pbar,              &
    & nC,nF,                  &
    & vCpn,                   &
    & tStoikioCpn,            &
    & tSimplex,iPosV,iZRov)
    ! & vFasIsPresent)

    call Simplex_StoreResult( &
    & iCount,                       &
    & TdgK,Pbar,                    &
    & size(vCpn),                   &
    & vCpn(:)%Mole,                 &
    & IPOSV,tSimplex(:,0))

    ! ! save results in table tSimplexResult
    ! call Simplex_StoreResult( &
    ! & iCount,TdgK,Pbar, &
    ! & IPOSV,tSimplex,   &
    ! & vFasIsPresent,tSimplexResult)

    if(iDebug>2) call Simplex_Print(vFas,size(vCpn),tSimplex(:,0),IPOSV)

  case(1)

    write(fTrc,'(A)') "Unbounded Objective Function"
    if(iDebug>0) &
    & print *,"SIMPLEX, Error: Unbounded Objective Function"

  case(-1)

    write(fTrc,'(A)') "No Solutions Satisfy Constraints Given"
    if(iDebug>0) &
    & print *,"SIMPLEX, Error: No Solutions Satisfy Constraints Given"

  case(-2)

    write(fTrc,'(A)') "BAD INPUT TABLEAU IN SIMPLX"
    if(iDebug>0) &
    & print *,"SIMPLEX, Error: BAD INPUT TABLEAU IN SIMPLX"

  end select

end subroutine Simplex_Run

subroutine Simplex_Print(vFas,nC,tSimplex_0,IPOSV)
  use M_T_Phase,     only: T_Phase
  !
  type(T_Phase),intent(in):: vFas(:)
  integer,      intent(in):: nC
  real(dp),     intent(in):: tSimplex_0(:)
  integer,      intent(in):: IPOSV(:)
  !
  integer:: i
  !
  do i=1,nC
    print '(G15.6,2X,A)',tSimplex_0(i),trim(vFas(IPOSV(i))%NamFs)
  end do
  print *,""
  !
end subroutine Simplex_Print

!-----------------------------------------------------------------------
!-- write tabulated results for all phases found stable at any step ----
!-----------------------------------------------------------------------
subroutine WriteSysComp(DimPath,vCpn)
  use M_Files,      only: DirOut
  use M_Dtb_Const,  only: T_CK
  use M_IoTools,    only: GetUnit
  use M_T_MixModel, only: MaxPole
  use M_T_Component,only: T_Component
  !
  use M_Global_Vars,only: vFas,vSpc,tFormula
  use M_Global_Vars,only: vMixModel,vDiscretModel,vDiscretParam
  !
  integer,          intent(in):: DimPath
  type(T_Component),intent(in):: vCpn(:)
  !
  integer :: I,J,K
  integer :: iPath,iFs,iCnt
  integer :: nC,nF
  integer :: iDis,iDisMod
  integer :: F
  logical,allocatable:: vDisIsPresent(:)
  integer,allocatable:: tDisParam(:,:,:)
  !type(T_Species):: Spc
  !
  nC=size(vCpn)
  nF=size(vFas)
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_moles.restab")
  !
  if(size(vDiscretModel)>0) then
    allocate(vDisIsPresent(size(vDiscretModel)))
    vDisIsPresent(:)= .false.
    do iFs=1,nF
     if(vFasIsPresent(iFs)) then
        !Spc= vSpc(vFas(iFs)%iSpc)
        iDis= vSpc(vFas(iFs)%iSpc)%iDiscret
        if(iDis/=0) &
        & vDisIsPresent(vDiscretParam(iDis)%iModel)= .true.
      end if
    end do
    !
    if(count(vDisIsPresent(:))>0) then
      ! N=0
      ! do J=1,size(vDiscretModel)
      !   if(vDisIsPresent(J)) &
      !   & N= N +vMixModel(vDiscretModel(J)%iMix)%NPole
      ! end do
      !print *,"dimension N=", N   ;  pause
      allocate(tDisParam(DimPath,3*size(vDiscretModel),3)) !MaxPole))
      tDisParam(:,:,:)= 0
    end if
  end if
  !
  !-----------------------------------------------------------title line
  write(F,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
  do iFs=1,nC
    write(F,'(A,A1)',advance="no") trim(vCpn(iFs)%NamCp)//"_cpn", T_
  end do
  ! do iFs=1,nC
  !   write(F,'(A,A1)',advance="no") trim(vCpn(iFs)%NamCp)//"_frac", T_
  ! end do
  do iFs=1,nF
    if(vFasIsPresent(iFs)) &
    & write(F,'(A,A1)',advance="no") trim(vFas(iFs)%NamFs), T_
  end do
  write(F,*)
  !----------------------------------------------------------/title line
  !
  iCnt=0
  do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then

      iCnt=iCnt+1
      write(F,'(I3,A1)',     advance="no") iCnt,T_


      write(F,'(2(F12.3,A1))',advance="no") &
      & tSimplexResult(nC+nF+1,iPath)-T_CK, T_, &
      & tSimplexResult(nC+nF+2,iPath),      T_

      !--- mole numbers of components
      do iFs=1,nC
        write(F,'(G15.6,A1)',advance="no") tSimplexResult(iFs,iPath),T_
      end do
      !---/mole numbers of components

      !~ Sum_=SUM(tSimplexResult(1:nC,iPath))
      !~ do iFs=1,nC
        !~ write(F,'(G15.6,A1)',advance="no") tSimplexResult(iFs,iPath)/Sum_,T_
      !~ end do

      do iFs=1,nF

        if(vFasIsPresent(iFs)) then

          write(F,'(G15.3,A1)',advance="no") &
          ! nr'moles phase          *nr'oxygen in formula
          & tSimplexResult(nC+iFs,iPath)*tFormula(iFs,1),T_

          iDis= vSpc(vFas(iFs)%iSpc)%iDiscret
          if(iDis/=0) then
            if(tSimplexResult(nC+iFs,iPath)>Zero) then
              iDisMod= vDiscretParam(iDis)%iModel
              !
              I= vDiscretParam(iDis)%I
              J= vDiscretParam(iDis)%J
              K= vDiscretParam(iDis)%K
              !
              if(I>=J .and. I>=K)      then  ;  iDisMod= iDisMod*3 -2
              elseif(J>=I .and. J>=K)  then  ;  iDisMod= iDisMod*3 -1
              else                           ;  iDisMod= iDisMod*3
              end if
              tDisParam(iCnt,iDisMod,1)= I
              tDisParam(iCnt,iDisMod,2)= J
              tDisParam(iCnt,iDisMod,3)= K
              !
              !write(15,'(3(I3,1X))',advance="NO") I,J,K
            end if
          end if

        end if

      end do
      write(F,*)

    end if

  end do
  !
  write(F,'(A1)') "_"
  !
  close(F)
  !
  if(allocated(tDisParam)) then
    call GetUnit(F)
    open(F,file=trim(DirOut)//"_mixtures.restab")
    !
    do I=1,iCnt
      do J=1,size(tDisParam,2)
        do K=1,3
          write(F,'(I3,1X)',advance="NO") tDisParam(I,J,K)
        end do
      end do
      write(F,*)
    end do
    deallocate(tDisParam)
    close(F)
  end if
  !
  if(allocated(vDisIsPresent)) deallocate(vDisIsPresent)
  !
  write(*,'(A)') "Results in ",trim(DirOut)//"_moles.restab"
  !
end subroutine WriteSysComp

subroutine Simplex_StoreResult( &
& iCount,                       &
& TdgK,Pbar,                    &
& nC,                           &
& vMolCpn,                      &
& IPOSV,tSimplex_0)
!--
!-- save results in table tSimplexResult
!--
  !! use M_Simplex_Vars,only: IPOSV,tSimplex
  !
  integer, intent(in):: iCount
  real(dp),intent(in):: TdgK,Pbar
  integer, intent(in):: nC
  real(dp),intent(in):: vMolCpn(1:nC)
  integer, intent(in):: IPOSV(1:nC)
  real(dp),intent(in):: tSimplex_0(0:nC+1)
  !
  integer::iC,iFs,N
  !
  N= size(tSimplexResult,1)
  !
  tSimplexResult(1:nC,iCount)= vMolCpn(1:nC)
  !
  do iC=1,nC
    iFs= IPOSV(iC)
    vFasIsPresent(iFs)= .true.
    tSimplexResult(nC+iFs,iCount)= tSimplex_0(iC) !! tSimplex(iC,0)
  end do
  !
  tSimplexResult(N-1,iCount)= TdgK
  tSimplexResult(N,  iCount)= Pbar
  !
end subroutine Simplex_StoreResult

subroutine Simplex_WriteResult( &
& TdgK,Pbar,                    &
& nCpn,nFs,                     &
& vCpn,                         &
& tStoikioCpn,                  &
& tSimplex,                     &
& IZROV,IPOSV)
!--
!-- write system, results, ...
!-- to files xxx_check1.log and xxx_check2.log
!--
  use M_IoTools
  use M_Dtb_Const,  only: T_CK
  use M_Files,      only: DirOut
  use M_T_Component,only: T_Component
  !
  use M_Global_Vars, only: vFas
  use M_Simplex_Build,only: fSpl1,fSpl2
  !
  real(dp),          intent(in):: TdgK,Pbar
  integer,           intent(in):: nCpn,nFs
  type(T_Component), intent(in):: vCpn(nCpn)
  real(dp),          intent(in):: tStoikioCpn(nFs,nCpn)
  real(dp),          intent(in):: tSimplex(0:nCpn+1,0:nFs)
  integer,           intent(in):: IZROV(1:nFs)
  integer,           intent(in):: IPOSV(1:nCpn)
  !
  integer,parameter:: maxOut1=100
  integer :: iCpn,iFas,I
  real(dp):: X
  ! character(len=30):: sFormat
  !
  ! nCpn= size(vCpn)
  ! nFs=  size(vFas)
  !
  if(fSpl2==0) then
    call GetUnit(fSpl2)
    open(fSpl2,file=trim(DirOut)//"_check2.log")
  end if
  !
  if(nFs<=maxOut1) then
    !
    ! file xxx_check2.log is written
    ! only when number of phases is not too big !!
    !
    if(fSpl1==0) then
      call GetUnit(fSpl1)
      open(fSpl1,file=trim(DirOut)//"_check1.log")
      write(fSpl1,'(4(A,/),/)') &
      & "line1= list of unstable phases", &
      & "line2= values for unstable phases", &
      & "lines 1..nC= stable phase name, amount, ...", &
      & "... stoikio of formation of unstable phase from stable assemblage"
    end if

    ! ! first line=
    ! ! list of unstable phases
    ! write(fSpl1,'(3(A,A1))',advance="no") "_",T_,"_",T_,"_",T_
    ! do iFas=1,nFs
    !   if(izrov(iFas)<=nFs) &
    !   & write(fSpl1,'(A,A1)',advance="no") trim(vFas(IZROV(iFas))%NamFs),T_
    ! end do
    ! write(fSpl1,*)
    !
    write(fSpl1,'(2(A,A1))',advance="no")   "_",T_,"_",T_
    write(fSpl1,'(G15.6,A1)',advance="no") tSimplex(0,0),T_

    ! second line=
    ! values for unstable phases
    do iFas=1,nFs
      if(izrov(iFas)<=nFs) &
      & write(fSpl1,'(G15.6,A1)',advance="no") tSimplex(0,iFas),T_
    end do
    write(fSpl1,*)
    !
    do iCpn=1,nCpn
      ! names of stable phases
      write(fSpl1,'(A,A1)',  advance="no") trim(vFas(IPOSV(iCpn))%NamFs),T_
      ! amount of stable phase
      write(fSpl1,'(G15.6,A1)', advance="no") tSimplex(iCpn,0),T_
      ! value for stable phase
      write(fSpl1,'(G15.6,A1)', advance="no") tSimplex(0,IPOSV(iCpn)),T_
      ! stoikio of formation of unstable phase from stable assemblage
      do iFas=1,nFs
        if (izrov(iFas)<=nFs) &
        & write(fSpl1,'(G15.6,A1)',advance="no") tSimplex(iCpn,iFas),T_
      end do !iFas
      write(fSpl1,*)
    end do !iCpn
    write(fSpl1,*)
    !
  end if
  !
  write(fSpl2,'(A)') "!"
  write(fSpl2,'(A,2G15.3)') "TdgC/PBar",TdgK-T_CK,Pbar
  write(fSpl2,'(A)') "!"
  !
  !----------------------------------------------------Write Composition
  write(fSpl2,'(A,A1)',advance="no") "vNamCpn",T_
  do iCpn=1,nCpn
    write(fSpl2,'(A,A1)',  advance="no") trim(vCpn(iCpn)%NamCp), T_
  end do
  write(fSpl2,*)
  !
  write(fSpl2,'(A,A1)',advance="no") "vMolCpn",T_
  do iCpn=1,nCpn
    write(fSpl2,'(G15.6,A1)',advance="no") vCpn(iCpn)%Mole, T_
  end do
  write(fSpl2,*)
  !
  write(fSpl2,'(A,A1)',  advance="no") "_",T_
  do iCpn=1,nCpn
    write(fSpl2,'(A,A1)',  advance="no") trim(vCpn(iCpn)%NamCp), T_
  end do
  write(fSpl2,*)
  !---------------------------------------------------/Write Composition

  ! !--------------------------------------------------Write Composition
  ! write(sFormat,'(A,I3,A)') '((A15,A1),',nCpn,'(A15,A1))'
  ! write(fSpl2,sFormat) "vNamCpn",T_,(vCpn(iCpn)%NamCp, T_),iCpn=1,nCpn)
 
  ! write(sFormat,'(A,I3,A)') '((A15,A1),',nCpn,'(G15.6,A1))'
  ! write(fSpl2,sFormat) "vMolCpn",T_,((vCpn(iCpn)%Mole,  T_),iCpn=1,nCpn)
 
  ! write(sFormat,'(A,I3,A)') '((A15,A1),',nCpn,'(A15,A1))'
  ! write(fSpl2,sFormat) "_",T_,      ((vCpn(iCpn)%NamCp,  T_),iCpn=1,nCpn)
  ! !-------------------------------------------------/Write Composition

  !---------------------------------------------------Check Mass Balance
  do iCpn=1,nCpn

    iFas=IPOSV(iCpn)
    write(fSpl2,'(A,A1)',  advance="no") trim(vFas(iFas)%NamFs),T_
    !
    vFasIsPresent(iFas)=.true.
    !
    do I=1,nCpn
      if( tStoikioCpn(iFas,I) /=Zero ) then
        write(fSpl2,'(G15.6,A1)',advance="no") &
        & tSimplex(iCpn,0)*tStoikioCpn(iFas,I), T_
      else
        write(fSpl2,'(A,A1)',advance="no") "0",T_
      end if
    end do
    write(fSpl2,*)

  end do
  !
  write(fSpl2,'(A)') "___"
  !
  write(fSpl2,'(A15,A1)',advance="no") "BALANCE",T_
  do I=1,nCpn
    !X=Zero
    !do iCpn=1,nCpn
    !  X=X + tSimplex(iCpn,0)*tStoikioCpn(IPOSV(iCpn),I)
    !end do
    X=SUM(tSimplex(1:nCpn,0)*tStoikioCpn(IPOSV(1:nCpn),I))
    write(fSpl2,'(G15.6,A1)',advance="no") X, T_
  end do
  write(fSpl2,*)
  !--------------------------------------------------/Check Mass Balance
  !
  return
end subroutine Simplex_WriteResult

end module M_Simplex_Path
