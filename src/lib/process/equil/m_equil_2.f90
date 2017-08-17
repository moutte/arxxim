module M_Equil_2
!--
!-- derived from M_Equil_3 --
!--
!-- calculate equilibrium speciation --
!-- considering saturation with minerals or gases --
!-- this version considers pure phases and mixtures --
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  implicit none

  private

  public:: Equil_Eq2

  real(dp),parameter:: FasMinim= 1.D-6
  ! FasMinim= mimim' amount of phase
  real(dp),parameter:: AffinIota= 1.D-3
  ! AffinIota= tolerance for equilibrium condition
  ! ( equilibrium ifF abs(Affinity)<AffinIota )
  real(dp),parameter:: MixMinim_TolX= 1.D-3
  ! for mixture minimisation (non ideal case)

  integer, allocatable:: tIPole(:,:)

  integer,parameter:: MaxIter= 30
  character(len=30):: sFMT

  character:: cTest_

contains

subroutine Equil_Eq2(NoMixture,iErr,cTest)
!--
!-- this version considers pure phases and mixtures
!-- (NoMixture - pure phase only )
!--
  use M_IOTools,      only: GetUnit
  use M_Files,        only: DirOut
  !
  use M_Numeric_Const,only: Ln10
  use M_Numeric_Tools,only: iMinLoc_R
  !
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: OutStrVec
  !
  use M_T_MixModel,   only: MaxPole
  use M_T_MixPhase,   only: T_MixPhase
  use M_T_Phase,      only: T_Phase
  !
  use M_Equil_Vars,   only: T_EquPhase
  use M_Equil_Solve,  only: Equil_Solve
  !
  use M_Global_Vars,  only: vSpc,vMixModel
  use M_Global_Vars,  only: vFas
  !
  use M_Basis_Vars,   only: tNuFas,vOrdPr
  use M_Equil_Vars,   only: cEquMode
  use M_Equil_Vars,   only: vYesList
  use M_Equil_Vars,   only: vDeltaG_Fs,vAffScale,vLnAct,vFasMole
  use M_Equil_Vars,   only: vEquFas,nEquFas
  use M_Equil_Vars,   only: vDeltaG_Eq,tNuEq
  !---------------------------------------------------------------------
  logical,intent(in) :: NoMixture
  integer,intent(out):: iErr
  character,intent(in),optional:: cTest
  !---------------------------------------------------------------------
  integer :: nCp,nFs,nPur,nSol
  integer :: iFs,iMix,iPur
  integer :: FF
  integer :: I,J,K,C
  integer :: nMix,nP !,nEquMix
  integer :: iPrecip,iElimin
  integer :: nIts,iDo0,iDo1
  ! real(dp):: X
  ! logical :: MixIsStable
  !
  integer, allocatable:: vIPole(:)
  real(dp),allocatable:: vNuFas(:)
  !
  real(dp),allocatable:: vFasAff(:)
  real(dp),allocatable:: vFasAffPur(:)
  !
  type(T_EquPhase),allocatable:: vEquFas0(:)
  type(T_MixPhase),allocatable:: vMixFas0(:)
  type(T_Phase),   allocatable:: vFas0(:) !,vFasPur(:)
  type(T_EquPhase):: EquFasNew
  !---------------------------------------------------------------------
  !! type:: T_EquPhase
  !!   character(len=23):: NamEq
  !!   integer :: iPur
  !!   integer :: iMix
  !!   real(dp):: vXPole(1:MaxPole)
  !!   !real(dp):: Grt
  !!   real(dp):: Mole
  !! end type T_EquPhase

  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Eq2_Mix"

  if(present(cTest)) then  ;  cTest_= cTest
  else                     ;  cTest_= "1"
  end if

  FF=0
  if(iDebug>1) then
    call GetUnit(FF)
    open(FF,file=trim(DirOut)//"_equil.log")
  end if

  cEquMode="EQ2"
  !-> flag in Equil_Solve, Equil_Residual, Equil_Jacobian

  nCp= size(vOrdPr)
  nPur= count(vFas(:)%iSpc/=0)
  nMix= size(vMixModel)
  nSol= 1 !! size(vSolFas)
  nFs= nPur +nMix !+ nSol

  !! DEBUGG !!
  if(iDebug>1) then
    write(6,'(A)') "==Equil_Eq2=="
    print *,"nPur=",nPur
    do J=1,size(vFas)
      !write(6,'(A)')  trim(vFas(J)%NamFs)
      print *,trim(vFas(J)%NamFs),vYesList(J)
    end do
  end if
  !pause

  if(NoMixture) nMix= 0

  allocate(vFasAff(nFs+nCp))
  allocate(vFasAffPur(nPur))
  allocate(vEquFas0(nCp))
  
  allocate(vFas0(nFs))
  vFas0(1:nPur)= vFas(1:nPur)
  
  !print *,"nSol, nPur, size(vFas), size(vFas0)=", &
  !&        nSol, nPur, size(vFas), size(vFas0)  ;  pause

  !if(nSol>0) vFas0(nPur+1)= vFas(nPur+1)
  
  !-------------------------------------------------------------mixtures
  if(nMix>0) then

    do K=1,nMix
      vFas0(nPur+nSol+K)%NamFs= vMixModel(K)%Name
    end do

    ! for each mixing model,
    ! we shall compute ONE composition of minimal G
    allocate(vMixFas0(nMix))
    do K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    end do

    allocate(tIPole(nMix,MaxPole))
    tIPole(:,:)= 0
    call MixModel_Init(nPur,vFas0,vMixModel)

  end if
  !------------------------------------------------------------/mixtures

  I=0
  do iPur=1,nPur
    if(vFas0(iPur)%MolFs>FasMinim) then
      I=I+1
      vEquFas(I)%iPur= iPur
      vEquFas(I)%iMix= 0
      vEquFas(I)%NamEq= trim(vFas0(iPur)%NamFs)
      vEquFas(I)%MolFs= vFas0(iPur)%MolFs
    end if
  end do
  nEquFas= I

  !------------------------------------------loop on mineral saturation,
  !----------------------------exit when all minerals have vPurYes false
  iDo0= 0
  Do0: do
    iDo0= iDo0 +1

    if(iDo0> MaxIter) then
      if(iDebug>2) then
        write(6,*) "iDo0>MaxIter in Equil_2"
        call Pause_
      end if
      iErr= -7
      exit Do0
    end if

    iDo1= 0
    !--------------------------solve the system until no negative phases
    Do1: do
      iDo1= iDo1 +1

      if(iDebug>2) then
        write(6,'(A)') "<-- vEquFas --"
        do I=1,nEquFas
          call EquFas_Trace(6,vEquFas(I))
        end do
        write(6,'(A)') "</- vEquFas --"
      end if

      if(nEquFas>0) then
        call EquFas_Clean
        call EquFas_Alloc(nCp,nEquFas)
        call EquFas_Update(nCp,nEquFas,vEquFas)
      end if

      !-------------------------------------------------solve the system
      call Equil_Solve(nIts,iErr)
      !----------------------------------------------------------------/

      if(nEquFas==0) exit Do1

      !---------------------------------------exclude "near zero" phases
      vEquFas0(:)= vEquFas(:)
      J= 0
      do I=1,nEquFas
        if(abs(vEquFas0(I)%MolFs)>FasMinim) then
          J=J+1
          vEquFas(J)= vEquFas0(I)
        end if
      end do
      nEquFas= J
      !--------------------------------------/exclude "near zero" phases
      
      !----------------------------check for phases with negative amount
      if(nEquFas>0) then
        !
        !----------------------------------------------------------trace
        if(iDebug>2) then
          write(6,'(A)') "<-- CURRENT --"
          do I=1,nEquFas
            call EquFas_Trace(6,vEquFas(I))
          end do
          write(6,'(A)') "</- CURRENT --"
        end if
        !---------------------------------------------------------/trace
        !
        !-------------------------remove most "negative" phase and cycle
        if(MINVAL(vEquFas(1:nEquFas)%MolFs) < -FasMinim) then

          iElimin= iMinLoc_R(vEquFas(1:nEquFas)%MolFs)

          if(iDebug>2) &
          & write(6,'(16X,2A)') "NEG= ", trim(vEquFas(iElimin)%NamEq)

          call EquFas_Remove(iElimin,nEquFas)

          cycle Do1 

        end if
        !------------------------/remove most "negative" phase and cycle
        !
      end if
      !
      exit Do1
      !---------------------------/check for phases with negative amount
      !
    end do Do1
    !-------------------------/solve the system until no negative phases

    if(iErr<0) exit Do0 != error in Newton -> exit

    vFasAff= Zero
    vFasAffPur= Zero

    do iPur=1,nPur

      if(vYesList(iPur)) &
      & vFasAffPur(iPur)= &
      & (  vDeltaG_Fs(iPur) &
      &  - dot_product(tNuFas(iPur,1:nCp),vLnAct(vOrdPr(1:nCp))) )

      ! scaling ala EQ3/6
      if(vYesList(iPur)) &
      & vFasAff(iPur)= vFasAffPur(iPur) /sum(abs(tNuFas(iPur,:))) !vAffScale(iPur)

    end do

    !-----------------------------------------------------------mixtures
    if(nMix>0) then

      do iMix=1,nMix

        ! for each mixing model,
        ! compute the composition that minimize
        ! the affinity of formation from end-members
        ! at current affinity w.r.t. the aqu'phase

        !vIPole(:)= tIPole(iMix,:)
        call Mixture_Minimize( &
        & vFasAffPur(:),tIPole(iMix,:),vMixModel(iMix), &
        & vMixFas0(iMix))

        vFasAff(nPur+iMix)= vMixFas0(iMix)%Grt

        !--------------------------------------------------------scaling
        nP= vMixModel(iMix)%NPole
        allocate(vNuFas(nCp))
        allocate(vIPole(nP))
        vIPole(1:nP)= tIPole(iMix,1:nP) !vMixModel(iMix)%vIPole(1:nP)
        do C=1,nCp
          vNuFas(C)= Zero
          do K=1,nP
            vNuFas(C)= vNuFas(C) &
            &        + vMixFas0(iMix)%vXPole(K) *tNuFas(vIPole(K),C)
          end do
          !! ! vNuFas(C)= sum( vMixFas0(iMix)%vXPole(1:nP) &
          !! ! &              * tNuFas(vIPole(1:nP),C) )
        end do
        deallocate(vIPole)
        vFasAff(nPur+iMix)= vMixFas0(iMix)%Grt /sum(abs(vNuFas(:)))
        deallocate(vNuFas)
        !-------------------------------------------------------/scaling

      end do

    end if
    !----------------------------------------------------------/mixtures

    !--------------------------------------------------------------trace
    if(iDebug>2) then
      write(6,'(A)') "<-- OVERSAT --"
      ! print *,size(vFasAff),nPur  ;  pause
      do I=1,nPur
        if (vYesList(I)) &
        ! .and. vFasAff(I)<AffinIota) &
        & write(6,'(A,G15.6,1X,A15,1X)') &
        & "Aff=",vFasAff(I),vFas0(I)%NamFs
      end do
      do I=1,nMix
        if (vFasAff(nPur+I)<AffinIota) then
          write(6,'(A,G15.6)',advance="NO") "Aff=",vFasAff(nPur+I)
          write(6,'(A)',advance="NO") " X(:)="
          do J=1,vMixModel(I)%nPole
            write(6,'(G15.6,1X)',advance="NO") vMixFas0(I)%vXPole(J)
          end do
          write(6,'(A)') vMixModel(I)%Name
        end if
      end do
      write(6,'(A)') "</- OVERSAT --"
      call Pause_
    end if
    !-------------------------------------------------------------/trace

    !---------------------add phase with highest supersaturation, if any
    !---------if no phase found with affinity below -AffinIota, exit Do0
    if(MINVAL(vFasAff)<-AffinIota) then

      iPrecip= iMinLoc_R(vFasAff)
      !-> the saturated phase with lowest vFasAff
      !! print *,"vFas(iPp)%NamFs= ",vFas(iPp)%NamFs  ;  pause
      !! vNewFas(iPrecip)= .true.

    else

      iPrecip= 0 ! there is no saturated phase,
      exit Do0   ! no phase found with affinity below -AffinIota -> exit

    end if

    !-------------------------------------------iPrecip>0, add new phase
    
    if(iPrecip<=nPur) then
    !--------------------------------------------------new phase is pure

      EquFasNew%iPur= iPrecip
      EquFasNew%iMix= 0
      EquFasNew%NamEq= trim(vFas0(iPrecip)%NamFs)
      EquFasNew%MolFs= Zero

    else
    !-----------------------------------------------new phase is mixture

      iMix= iPrecip -nPur

      EquFasNew%iPur= 0
      EquFasNew%iMix= iMix
      EquFasNew%NamEq= trim(vMixModel(iMix)%Name)
      nP= vMixModel(iMix)%NPole
      EquFasNew%nPole= nP
      EquFasNew%vIPole(1:nP)= tIPole(iMix,1:nP)

      EquFasNew%vXPole(1:nP)= vMixFas0(iMix)%vXPole(1:nP)
      EquFasNew%vLnAct(1:nP)= vMixFas0(iMix)%vLnAct(1:nP)
      EquFasNew%MolFs= Zero

    end if

    call CheckAssemblage(FF,nCp,nEquFas,EquFasNew,vEquFas,iElimin)

    !----------------------if necessary, remove one of the active phases
    if(iElimin>0) then

      call EquFas_Remove(iElimin,nEquFas)

      if(iDebug>2) &
      & write(6,'(16X,2A)') "OUT= ",trim(vEquFas(iElimin)%NamEq)

    end if
    !-----------------------------------/remove one of the active phases

    vEquFas(nEquFas+1)= EquFasNew
    nEquFas= nEquFas +1

    if(iDebug>2) &
    & write(6,'(16X,2A)') "ADD= ",trim(EquFasNew%NamEq)
    !-----------------------------------------------------/add new phase
    
    !--------------------------------------------------------------trace
    if(iDebug>3) then
      do I=1,nEquFas
        iFs= vEquFas(I)%iPur
        if(iFs>0) then
          write(6,'(A,I4,1X,A15,2(A,G15.6))') &
          & "Do0",iDo0, &
          & vFas0(iFs)%NamFs, &
          & "/ lQsK=", -vFasAff(iFs)/Ln10 ,&
          & "/ Mole=",vEquFas(I)%MolFs
        end if
      end do
      !! if(iDebug>3)
      call Pause_
    end if
    !-------------------------------------------------------------/trace
    !
    !-------------------/ add phase with highest supersaturation, if any
    !
    if(nEquFas==0) exit Do0
    !
    if(iDebug>2) print *,"============================================="
    !! if(iDebug>2) pause

  end do Do0
  !------------------------------------------/loop on mineral saturation

  !! if (allocated(vFasMole)) deallocate(vFasMole)
  !! allocate(vFasMole(nPur+nMix))
  !! vFasMole(:)= Zero
  !! do I=1,nEquFas
  !!   if(vEquFas(I)%iPur>0) vFasMole(vEquFas(I)%iPur)= vEquFas(I)%Mole
  !!   if(vEquFas(I)%iMix>0) vFasMole(vEquFas(I)%iMix+nPur)= vEquFas(I)%Mole
  !! end do

  if(nMix>0) call Check_vXMean

  call EquFas_Clean

  deallocate(vFas0)
  deallocate(vFasAff)
  deallocate(vFasAffPur)
  deallocate(vEquFas0)

  if(nMix>0) then
    deallocate(vMixFas0)
    deallocate(tIPole)
  end if

  if(FF>0) close(FF)

  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Eq2_Mix"
  !
contains

subroutine Check_vXMean
!--
!-- for each mixing model with nFas>1
!-- compute the free energy of the mixture
!-- with the mean composition
!--
  use M_System_Vars,only: TdgK,Pbar
  !
  !--- T_TabPhase: same type as used for GEM routine:
  ! for a given mixing model,
  ! there can be, in the stable assemblage,
  ! up to NPole phases of different compositions.
  ! T_TabPhase is design to contain the different compositions
  ! that may be active for a same mixing model.
  ! nFas is the number of active phases
  type:: T_TabPhase
    character(len=23):: NamTbl
    integer :: nFas                        ! -> number of phases
    real(dp):: tXPole(1:Maxpole,1:MaxPole) ! table of phase compositions
    real(dp):: vMole(1:MaxPole)            ! mole nr of each phase
  end type T_TabPhase
  type(T_TabPhase),allocatable:: vTabPhase(:)
  !---/(same as for GEM routine)
  !
  type(T_EquPhase):: EquFas
  integer :: I,J,K,P,nP,N,nF
  integer :: nTabFas
  real(dp):: G
  logical :: MustUpdate
  real(dp),allocatable:: vX(:)
  real(dp),allocatable:: vG0(:)

  MustUpdate= .false.
  allocate(vTabPhase(nMix))
  vTabPhase(:)%nFas= 0
  nTabFas= 0

  do I=1,nEquFas
    EquFas= vEquFas(I)
    if(EquFas%iMix>0) then
      J= EquFas%iMix
      vTabPhase(J)%NamTbl= trim(EquFas%NamEq)
      vTabPhase(J)%nFas= vTabPhase(J)%nFas +1
      vTabPhase(J)%vMole(vTabPhase(J)%nFas)= EquFas%MolFs
      nP= vMixModel(J)%NPole
      do P=1,nP
        vTabPhase(J)%tXPole(vTabPhase(J)%nFas,P)= EquFas%vXPole(P)
      end do
    end if
  end do

  if(count(vTabPhase(:)%nFas>0)==0) return

  do J=1,nMix

    nF= vTabPhase(J)%nFas
    if (nF>1) then

      nP= vMixModel(J)%NPole

      !-- compute the mean composition of the mixture
      allocate(vX(nP))
      vX(:)= Zero
      do K=1,nF
        vX(1:nP)= vX(1:nP) &
        &       + vTabPhase(J)%tXPole(K,1:nP) &
        &         *vTabPhase(J)%vMole(K)
      end do
      vX(:)= vX(:) /sum(vX(:))

      !-- compute its free energy,
      !-- using vFasAffPur(:),tIPole(J,:),vMixModel(J)
      allocate(vG0(nP))
      do P=1,nP
        vG0(P)= vFasAffPur(tIPole(J,P))
      end do
      !
      G= MixPhase_Grt( &
      & TdgK,Pbar, &
      & nP, &
      & vMixModel(J), &
      & vG0(1:nP), &
      & vX(1:nP))
      !
      print *,"G=",G
      ! if G is near zero,
      ! then there is only one phase, with the composition vX
      if(abs(G)<AffinIota *1.D1) then
        MustUpdate= .true.
        vTabPhase(J)%nFas= 1
        vTabPhase(J)%tXPole(1,1:nP)= vX(1:nP)
        vTabPhase(J)%vMole(1)= sum(vTabPhase(J)%vMole(1:nF))
      end if

      deallocate(vG0)
      deallocate(vX)

    end if

  end do

  if(MustUpdate) then

    K= 0
    vEquFas0(:)= vEquFas(:)
    vEquFas(:)%MolFs= Zero

    do I=1,nEquFas
      if(vEquFas0(I)%iMix==0) then
        K=K+1
        vEquFas(K)= vEquFas0(I)
      end if
    end do

    do J=1,nMix
      if(vTabPhase(J)%nFas>0) then
        !! ! sum( vTabPhase(J)%vMole(1:vTabPhase(J)%nFas) )
        !! J= vEquFas0(I)%iMix
        !! print *, "J=", J  ;  pause
        do N= 1,vTabPhase(J)%nFas
          K=K+1
          vEquFas(K)%NamEq= trim(vTabPhase(J)%NamTbl)
          vEquFas(K)%iPur= 0
          vEquFas(K)%iMix= J
          vEquFas(K)%vXPole(1:nP)= vTabPhase(J)%tXPole(N,1:nP)
          vEquFas(K)%MolFs= vTabPhase(J)%vMole(N)
        end do
      end if
    end do

    nEquFas= K

  end if

  deallocate(vTabPhase)

  return
end subroutine Check_vXMean

subroutine EquFas_Remove(N,M)
  integer,intent(in)   :: N
  integer,intent(inout):: M

  integer:: I,J

  vEquFas0(:)= vEquFas(:)
  J= 0
  do I=1,M
    if(I/=N) then
      J=J+1
      vEquFas(J)= vEquFas0(I)
    end if
  end do
  M= J

  return
end subroutine EquFas_Remove

end subroutine Equil_Eq2

subroutine EquFas_Alloc(nC,nEq)
  use M_Equil_Vars,only: vDeltaG_Eq,tAlfEq,tNuEq
  !
  integer,intent(in):: nC,nEq
  !
  allocate(vDeltaG_Eq(nEq))
  allocate(tAlfEq(nC,nEq))
  allocate(tNuEq(nEq,nC))
  !
  return
end subroutine EquFas_Alloc

subroutine EquFas_Update(nCp,nEquFas,vEquFas)
!--
!-- update vDeltaG_Eq, tAlfEq, tNuEq,
!-- which are used in Equ_Residual
!--
  use M_T_MixModel, only: MaxPole
  use M_Global_Vars,only: vMixModel
  use M_Equil_Vars, only: T_EquPhase
  use M_Basis_Vars, only: tAlfFs,tNuFas
  use M_Equil_Vars, only: vDeltaG_Eq,tAlfEq,tNuEq
  use M_Equil_Vars, only: vDeltaG_Fs
  !---------------------------------------------------------------------
  integer,intent(in):: nCp
  integer,intent(in):: nEquFas
  type(T_EquPhase),intent(in):: vEquFas(:)
  !---------------------------------------------------------------------
  integer:: iPur,iMix,I,C

  do I=1,nEquFas

    iPur= vEquFas(I)%iPur
    iMix= vEquFas(I)%iMix

    if(iPur>0) then
      vDeltaG_Eq(I)= vDeltaG_Fs(iPur)
      tAlfEq(:,I)= tAlfFs(:,iPur)
      tNuEq(I,:)= tNuFas(iPur,:)
    end if

    if(iMix>0) call EquFasMix_Update(I,nCp,vEquFas(I))

  end do

  if(iDebug>2) then
    ! write(sFMT,'(a,i3,a)') '(A,1X,',nCp,'(G12.3,1X))'
    write(75,'(A)') "===========================EquFas_Update=========="
    do I=1,nEquFas
      write(75,'(A,G15.6,1X,A)') &
      & "DeltaG_Eq= ",vDeltaG_Eq(I),trim(vEquFas(I)%NamEq)
      !write(75,sFMT) "tAlfEq=    ",(tAlfEq(C,I),C=1,nCp)
      !write(75,sFMT) "tNuEq=     ",(tNuEq(I,C),C=1,nCp)
      write(75,'(A,*(G12.3,1X))') "tAlfEq=    ",(tAlfEq(C,I),C=1,nCp)
      write(75,'(A,*(G12.3,1X))') "tNuEq=     ",(tNuEq(I,C),C=1,nCp)
    end do
  end if
  !
  return
end subroutine EquFas_Update

subroutine EquFasMix_Update(I,nCp,EquFas)
!--
!-- compute vDeltaG_Eq(I),tAlfEq(:,I),tNuEq(I,:)
!-- using EquFas%vXPole(:), and vDeltaG_Fs,tAlfFs,tNuFas
!--
  use M_Equil_Vars, only: T_EquPhase
  use M_Basis_Vars, only: tAlfFs,tNuFas  !IN
  use M_Equil_Vars, only: vDeltaG_Fs     !IN
  use M_Equil_Vars, only: vDeltaG_Eq,tAlfEq,tNuEq !OUT
  !---------------------------------------------------------------------
  integer,intent(in):: I
  integer,intent(in):: nCp
  type(T_EquPhase),intent(in):: EquFas
  !---------------------------------------------------------------------
  integer:: C,nP
  integer,allocatable:: vIPole(:) ! indexes of end-members in vFas

  nP= EquFas%NPole
  allocate(vIPole(nP))
  vIPole(1:nP)= EquFas%vIPole(1:nP)

  vDeltaG_Eq(I)= sum( EquFas%vXPole(1:nP) &
  &                 * vDeltaG_Fs(vIPole(1:nP)) ) &
  &            + sum( EquFas%vXPole(1:nP) &
  &                 * EquFas%vLnAct(1:nP) )

  do C=1,nCp
    tAlfEq(C,I)= sum( EquFas%vXPole(1:nP) &
    &               * tAlfFs(C,vIPole(1:nP)) )
    tNuEq(I,C)= sum( EquFas%vXPole(1:nP) &
    &              * tNuFas(vIPole(1:nP),C) )
  end do

  deallocate(vIPole)

  return
end subroutine EquFasMix_Update

subroutine EquFas_Trace(iFile,Fas)
  use M_Equil_Vars,only: T_EquPhase
  !
  integer,         intent(in):: iFile
  type(T_EquPhase),intent(in):: Fas
  !
  integer:: J
  !
  write(iFile,'(A,G15.6,1X)',advance="NO") "N=",Fas%MolFs
  if(Fas%iMix>0) then
    write(iFile,'(A)',advance="NO") " X(:)="
    do J=1,Fas%nPole
      write(iFile,'(G15.6,1X)',advance="NO") Fas%vXPole(J)
    end do
  end if
  write(iFile,'(A)') trim(Fas%NamEq)
  !
end subroutine EquFas_Trace

subroutine EquFas_Clean
  use M_Equil_Vars,only: vDeltaG_Eq,tAlfEq,tNuEq

  if(allocated(vDeltaG_Eq)) deallocate(vDeltaG_Eq)
  if(allocated(tAlfEq))     deallocate(tAlfEq)
  if(allocated(tNuEq))      deallocate(tNuEq)

  return
end subroutine EquFas_Clean

subroutine CheckAssemblage( & ! IN
& FF,                       & ! IN
& nCp,                      & ! IN
& nEquFas,                  & ! IN
& EquFasNew,                & ! IN, the new phase
& vEquFas,                  & ! IN, the current phase set
& iElimin)                    ! OUT
!--
!-- check whether new species, indexed vFas(iPp)%iPur
!-- is independent of current phase assemblage
!--
  use M_Numeric_Mat,  only: LU_BakSub
  use M_Numeric_Tools,only: iMinLoc_R,iFirstLoc,iMaxLoc_R
  use M_Basis_Tools,  only: Basis_FreeSet_Select
  use M_T_MixModel,   only: MaxPole
  use M_Equil_Vars,   only: T_EquPhase
  !
  use M_Basis_Vars,   only: tAlfPr,tAlfFs,vOrdPr
  !---------------------------------------------------------------------
  integer, intent(in)   :: FF
  integer, intent(in)   :: nCp
  integer, intent(in)   :: nEquFas
  type(T_EquPhase),intent(in):: EquFasNew
  type(T_EquPhase),intent(in):: vEquFas(:)
  integer, intent(out)  :: iElimin
  !---------------------------------------------------------------------
  integer :: iPur,iMix,I,J,C,nP
  logical :: IsNotFree
  logical :: bSingul
  integer :: vIPole(MaxPole)
  !
  ! integer :: iElimin1, iElimin2
  ! real(dp):: xMin, xMax
  !
  real(dp),allocatable:: tStoikCpn(:,:),tStoikCpn2(:,:)
  real(dp),allocatable:: tBase(:,:)
  integer, allocatable:: vIndex(:)
  integer, allocatable:: vIndx(:)
  real(dp),allocatable:: vY(:)
  real(dp),allocatable:: vTest(:),vTest2(:)
  !---------------------------------------------------------------------
  !! if(iDebug>2) write(6,'(A)') "< CheckAssemblage"
  !
  if(iDebug>2) write(FF,'(2A)') "NEW= ", EquFasNew%NamEq
  !
  iElimin= 0
  !
  allocate(vIndex(nCp))
  vIndex(:)=0
  do i=1,nEquFas+2
    vIndex(i)= i
  end do
  allocate(tStoikCpn(nCp,nEquFas+2))
  !
  !----------------------------------------------------compute tStoikCpn
  !----------------------= stoikio table of current set of stable phases
  tStoikCpn(:,1)= tAlfPr(:,1) ! water
  vIndex(1)= 1
  !
  do I=1,nEquFas

    iPur= vEquFas(I)%iPur
    iMix= vEquFas(I)%iMix

    if(iPur>0) tStoikCpn(:,1+I)= tAlfFs(:,iPur)

    if(iMix>0) then
      nP= vEquFas(I)%nPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      do C=1,nCp
        tStoikCpn(C,1+I)= sum( vEquFas(I)%vXPole(1:nP) &
        &                    * tAlfFs(C,vIPole(1:nP)) )
      end do
    end if

  end do

  iPur= EquFasNew%iPur
  iMix= EquFasNew%iMix

  if(iPur>0) tStoikCpn(:,nEquFas+2)= tAlfFs(:,EquFasNew%iPur)

  if(iMix>0) then
    nP= EquFasNew%nPole
    vIPole(1:nP)= EquFasNew%vIPole(1:nP)
    do C=1,nCp
      tStoikCpn(C,nEquFas+2)= sum( EquFasNew%vXPole(1:nP) &
      &                          * tAlfFs(C,vIPole(1:nP)) )
    end do
  end if
  !---------------------------------------------------/compute tStoikCpn

  !----------------------------------------------------------------trace
  if(iDebug>2) then
    do i=1,nEquFas
      do C=1,nCp
        write(71,'(F7.3,1X)',advance="no") tStoikCpn(C,i+1)
      end do
      write(71,'(A)') trim(vEquFas(i)%NamEq)
    end do
    do C=1,nCp
      write(71,'(F7.3,1X)',advance="no") tStoikCpn(C,nEquFas+2)
    end do
    write(71,'(A)') trim(EquFasNew%NamEq)
    write(71,'(A)') "=================================================="
  end if
  !---------------------------------------------------------------/trace

  !-- using Gaussian elimination (without permutations),
  !-- find a set of independent species (of stoikios tAlfPr)
  !-- independent also of those in tStoikioCpn
  call Basis_FreeSet_Select( & !
  & tAlfPr,      & !IN
  & vIndex,      & !OUT
  & IsNotFree,   & !OUT
  & tStoikioCpn= tStoikCpn)     !IN
  !
  deallocate(vIndex)

  if(IsNotFree) then
    !-- new phase not independent from current assemblage
    !-- -> build a basis without the new phase,
    !-- and compute stoichio of new phase versus current assemblage
    !
    allocate(tStoikCpn2(nCp,nEquFas+1))
    tStoikCpn2(:,1:nEquFas+1)= tStoikCpn(:,1:nEquFas+1)

    allocate(vIndex(nCp))
    vIndex(:)=0
    do i=1,nEquFas+1
      vIndex(i)= i
    end do

    !-- find a set of independent species (of stoikios tAlfPr)
    !-- independent also of those in tStoikioCpn2
    !-- vIndex(:) gives the indexes in tAlfPr of these independent species
    call Basis_FreeSet_Select( & !
    & tAlfPr,                  & !IN
    & vIndex,                  & !OUT
    & IsNotFree,               & !OUT
    & tStoikioCpn= tStoikCpn2)   !IN

    allocate(tBase(nCp,nCp))
    tBase(:,1:nEquFas+1)= tStoikCpn(:,1:nEquFas+1)
    do I=nEquFas+2,nCp
      tBase(:,I)= tAlfPr(:,vIndex(I))
    end do

    !--------------------------------------------------------------trace
    if(iDebug>2) then
      write(73,'(A)') "==CheckAssemblage=="
      do I=1,nCp
        do J=1,nCp
          write(73,'(G11.2,1X)',advance="NO") tBase(J,I)
        end do
        write(73,*)
      end do
      write(73,'(A)') "==**=="
    end if
    !-------------------------------------------------------------/trace

    allocate(vIndx(nCp))
    allocate(vY(nCp))

    !------------------------------------------LU decomposition of tBase
    call Compute_Transform(tBase,vIndx,bSingul)
    if(bSingul) call Stop_("SINGUL IN CheckAssemblage")

    iPur= EquFasNew%iPur
    iMix= EquFasNew%iMix
    if(iPur>0) vY(:)= tAlfFs(:,iPur)
    if(iMix>0) then
      nP= EquFasNew%nPole
      vIPole(1:nP)= EquFasNew%vIPole(1:nP)
      do C=1,nCp
        vY(C)= sum( EquFasNew%vXPole(1:nP) &
        &         * tAlfFs(C,vIPole(1:nP)) )
      end do
    end if
    call LU_BakSub(tBase,vIndx,vY)

    !-----------------------------------find which phase must be removed
    allocate(vTest(nEquFas))   ;  vTest= Zero
    allocate(vTest2(nEquFas))  ;  vTest2= 1.D9
    
    do I=1,nEquFas
      if(abs(vY(I+1))>1.D-9) vTest2(I)= vEquFas(I)%MolFs/vY(I+1)
      vTest(I)= abs(vY(I+1)) /vEquFas(I)%MolFs
    end do

    !--------------------------------------------------------------trace
    if(iDebug>2) then
      write(73,'(2A)') "EquFasNew=",trim(EquFasNew%NamEq)
      do I=1,nEquFas
        write(73,'(3(G15.6,1X),A)') &
        & vY(I+1),vEquFas(I)%MolFs,vTest(I),trim(vEquFas(I)%NamEq)
        !! write(73,'(2(G15.6,1X),A)') vY(I+1),vTest(I),trim(vEquFas(I)%NamEq)
      end do
    end if
    !-------------------------------------------------------------/trace

    !!X= maxval(vTest)
    select case(cTest_)
    case("2")  ;  iElimin= iMinLoc_R(vTest2)
    case("1")  ;  iElimin= iMaxLoc_R(vTest)
    end select

    deallocate(vTest)
    deallocate(vTest2)

    !----------------------------------/find which phase must be removed
    !
    deallocate(vIndx)
    deallocate(vY)
    deallocate(tBase)
    !
    deallocate(vIndex)
    deallocate(tStoikCpn2)
    !
  end if
  !
  deallocate(tStoikCpn)
  !
  return
end subroutine CheckAssemblage

subroutine MixModel_Init(nPur,vFas,vMixModel)
!--
!-- initialize tIPole(iMix,iP),
!-- - address of e-m iP of mixture iMix in vSpc or vFas
!--
  use M_T_Phase,    only: T_Phase
  use M_T_MixModel, only: T_MixModel

  integer,         intent(in):: nPur
  type(T_Phase),   intent(in):: vFas(:)
  type(T_MixModel),intent(in):: vMixModel(:)

  integer:: iModel,P,P1,P2,iPur,nP
  !
  do iModel=1,size(vMixModel)

    nP= vMixModel(iModel)%NPole

    ! find the indexes of end-members of mixture vMixFas(iMix)
    ! in the pure phase list vFas(1:nPur)
    ! -> indices in tIPole(iMix,:)
    !! print *,trim(vMixModel(iModel)%Name)
    do P=1,nP
      P1= vMixModel(iModel)%vIPole(P)
      !-> index in vSpc
      !-> find the pure phases that points to this species
      P2= 0
      do iPur=1,nPur
        if(vFas(iPur)%iSpc==P1) then
          P2= iPur
          !! vPurIsPole(iPur)= .true.
          exit
        end if
      end do
      tIPole(iModel,P)= P2
      !! print *,trim(vFas(P2)%NamFs)
      !
    end do

  end do
  !
  return
end subroutine MixModel_Init

subroutine Compute_Transform(tTransform,vIndx,Error)
  use M_Numeric_Mat,only: LU_Decomp
  !
  real(dp),intent(inout):: tTransform(:,:)
  integer, intent(out)  :: vIndx(:)
  logical, intent(out)  :: Error
  !
  real(dp):: D
  !
  ! the transformation matrix, tTransform, must be invertible !!
  call LU_Decomp(tTransform,vIndx,D,Error)
  !
end subroutine Compute_Transform

real(dp) function MixPhase_Grt(TdgK,Pbar,nP,MM,vMu0rt,vX)
  use M_T_MixModel,only: T_MixModel,MixModel_Activities
  !----------------------------------------------------------------inout
  real(dp),        intent(in):: TdgK,Pbar
  integer,         intent(in):: nP
  type(T_MixModel),intent(in):: MM        ! mixing model
  real(dp),        intent(in):: vMu0rt(:) !
  real(dp),        intent(in):: vX(:)     ! phase composition
  !-----------------------------------------------------------------vars
  real(dp):: vLGam(nP),vLIdeal(nP),vLnAct(nP)
  logical :: vLPole(nP)
  real(dp):: G
  integer :: i
  logical :: Ok
  character(len=30):: Msg
  !---------------------------------------------------------------------
  vLPole(:)= (vX(:)>Zero)
  !
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
  !
  return
end function MixPhase_Grt

subroutine Mixture_Minimize( &
& vAffPole,vIPole,MixModel, &
& MixFas)
!--
!-- for a given mixture model, MixModel,
!-- and given affinities of its end-members, vAffPole,
!-- compute the composition X that minimizes G: saved in MixFas%vXPole(:)
!-- and compute G at X: saved in MixFas%Grt
!--
  use M_Numeric_Const, only: Ln10
  use M_Safe_Functions,only: FSafe_Exp
  use M_System_Vars,   only: TdgK,Pbar
  !
  use M_T_Phase,   only: T_Phase
  use M_T_MixPhase,only: T_MixPhase
  use M_T_MixModel,only: T_MixModel,MaxPole,Mix_Molecular
  !
  use M_Optimsolver_Theriak
  use M_MixModel_Optim
  !
  ! use M_Global_Vars,only: vFas
  !---------------------------------------------------------------------
  real(dp),        intent(in)   :: vAffPole(:)
  integer,         intent(in)   :: vIPole(:)
  type(T_MixModel),intent(in)   :: MixModel
  type(T_MixPhase),intent(inout):: MixFas
  !---------------------------------------------------------------------
  real(dp):: vX(MaxPole),vMu(MaxPole),vLnAct(MaxPole)
  real(dp):: vXmin(MaxPole),vMuMin(MaxPole)
  real(dp):: G,Gmin
  integer :: nP,K,P
  integer :: FF
  integer :: Multi
  !
  real(dp):: TolX,DeltaInit
  integer :: its,nCallG
  logical :: OkConverge

  FF= 0 ! log file
  TolX= MixMinim_TolX
  DeltaInit= 0.05D0

  nP= MixModel%NPole
  Multi= MixModel%vMulti(1)
  !
  !! !---- compute affin' of all reactions between end-memb' and water --
  !! do P= 1,nP
  !!   ! index of pole P in pure phase list
  !!   ! K= tIPole(iModel,P)
  !!   K= MixModel%vIPole(P)
  !!   ! compute affinity of reaction (pole P <> water )
  !!   !vAffPole(P)= vFas(K)%Grt &
  !!   !&          - dot_product(tNuFas(K,:), vSpc(vOrdPr(:))%G0rt) &
  !!   !&          - dot_product(tNuFas(K,:), vLnAct(vOrdPr(:)))
  !!   vAffPole(P)= vDeltaGPole(K) &
  !!   &          - dot_product(tNuFas(K,:), vLnAct(vOrdPr(:)))
  !! end do
  !! !---
  !
  if( MixModel%Model==Mix_Molecular .and. MixModel%NMarg==0 ) then
    !---------------------------------------------------analytic minimum
    Multi= MixModel%vMulti(1)
    Multi= 1

    do P=1,nP
      ! K= MixModel%vIPole(P)
      ! vXmin(P)= FSafe_Exp(-vFasPur(MixModel%vIPole(P))%Grt /real(Multi))
      vXmin(P)= FSafe_Exp(-vAffPole(vIPole(P)) /real(Multi))
    end do
    vXmin(1:nP)=  vXmin(1:nP) /sum(vXmin(1:nP))
    vLnAct(1:nP)= Multi*log(vXmin(1:nP))
    vMuMin(1:nP)= vAffPole(vIPole(1:nP)) + vLnAct(1:nP)

    Gmin= sum(vXmin(1:nP) * vMuMin(1:nP))

    !--------------------------------------------------------------trace
    !if(iDebug>2) then
    !write(74,'(A)') "--Mixture_Minimize--"
    !  do P=1,nP
    !    write(74,'(A,2G15.6,1X,A)') &
    !    & "vAffPole(P),vXmin(P)= ",&
    !    & vAffPole(vIPole(P)), &
    !    & vXmin(P), &
    !    & trim(vFas(vIPole(P))%NamFs)
    !  end do
    !end if
    !-------------------------------------------------------------/trace
    !--------------------------------------------------/analytic minimum
  else
    !-----------------------------------------------numerical minimum(s)
    call MixModel_Optim_SetParams(TdgK,Pbar,MixModel)
    allocate(Mixmodel_Optim_vMu0rt(nP))
    allocate(Mixmodel_Optim_vLPole(nP))
    ! Mixmodel_Optim_vMu0rt(1:nP)= vFasPur(MixModel%vIPole(1:nP))%Grt
    Mixmodel_Optim_vMu0rt(1:nP)= vAffPole(vIPole(1:nP))
    Mixmodel_Optim_vLPole(1:nP)= .true.

    Gmin= 1.0D30

    do P=1,nP

      !--- start from compos'nP close to end-member P
      vX(P)= One - 1.0D-3
      do K=1,nP
        if(K/=P) vX(K)= 1.0D-3/real(nP-1)
      end do

      !! vX(1:nP)= vMixFas_Xpole_Init(I)%tXPole(P,1:nP)

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

      !! if(.not. OkConverge) then
        !! print *,"NoConvergence in Mixture_Minimize"
        !! cycle
      !! end if
      !vMixFas_Xpole_Init(I)%tXPole(P,1:nP)= vX(1:nP)

      if(G<Gmin) then
        Gmin= G
        vXmin(1:nP)= vX(1:nP)
        vMuMin(1:nP)= vMu(1:nP)
      end if

    end do
    !pause

    deallocate(Mixmodel_Optim_vMu0rt)
    deallocate(Mixmodel_Optim_vLPole)
    !----------------------------------------------/numerical minimum(s)
  end if

  MixFas%vLPole(:)= .true.
  ! write(76,'(G15.6)') sum(vXmin(:))
  vXmin(:)= vXmin(:) /sum(vXmin(:))
  MixFas%vXPole(:)= vXmin(:)
  MixFas%Grt= Gmin
  MixFas%vLnAct(1:nP)= vMuMin(1:nP) -vAffPole(vIPole(1:nP))

  !------------------------------------------------------------log files
  !! if(F1>0) then
  !!   write(F1,'(/,A,/)') "Mixture_Minimize,EndMember,Xi,Gi"
  !!   do I=1,size(vMixModel)
  !!     !if(vMixFas0(I)%Grt > 1.D-3) cycle
  !!     MixModel= vMixModel(I)
  !!     write(F1,'(2A)') "MODEL= ",MixModel%Name
  !!     do K=1,MixModel%NPole
  !!       write(F1,'(A,1X,G15.6)') &
  !!       & MixModel%vNamPole(K),        &
  !!       & vMixFas0(I)%vXPole(K)
  !!     end do
  !!     write(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
  !!     write(F1,*)
  !!   end do
  !!   write(F1,*)
  !! end if
  !! !
  !! if(F2>0) then
  !!   do I=1,size(vMixModel)
  !!   !if(vMixFas0(I)%Grt < Zero) then
  !!     MixModel= vMixModel(I)
  !!     !
  !!     write(F2,'(A,A1)',advance="NO") MixModel%Name,T_
  !!     do K=1,MixModel%NPole
  !!       write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%vXPole(K),T_
  !!     end do
  !!     !do K=1,MixModel%NPole
  !!     !  write(F2,'(G15.6,A1)',advance="NO") vFasPur(MixModel%vIPole(K))%Grt,T_
  !!     !end do
  !!     write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%Grt,T_
  !!     !
  !!   !end if
  !!   end do
  !!   write(F2,*)
  !! end if
  !-----------------------------------------------------------/log files
  !
end subroutine Mixture_Minimize

end module M_Equil_2
