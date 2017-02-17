module M_Equil_1
!--
!-- calculate equilibrium speciation --
!-- considering saturation with minerals or gases --
!-- this version considers pure phases and mixtures --
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  implicit none
  !
  private
  !
  public:: Equil_Eq1
  !
  real(dp),parameter:: MixMinim_TolX= 1.D-4
  != parameter for Mixture_Minimize
  !
  real(dp),parameter:: FasMinim= 1.D-9
  ! FasMinim= mimim' amount of phase
  !
  real(dp),parameter:: AffinIota= 1.D-6
  ! AffinIota= tolerance for equilibrium condition
  ! ( equilibrium == ABS(Affinity)<AffinIota )
  !
  real(dp),parameter:: dXi_Minim= 1.D-16
  real(dp),parameter:: dXi_Maxim= 1.D+12
  !
  integer, allocatable:: tIPole(:,:)
  character(len=30):: sFMT

contains

subroutine Equil_Eq1(NoMixture,iErr)
!--
!-------------------- this version considers pure phases and mixtures --
!--
  use M_IOTools,      only: GetUnit
  use M_Files,        only: DirOut
  use M_Numeric_Const,only: Ln10
  use M_Numeric_Tools,only: iMinLoc_R
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: OutStrVec
  use M_T_MixModel,   only: MaxPole
  use M_T_MixPhase,   only: T_MixPhase
  use M_T_Phase,      only: T_Phase
  !
  use M_Equil_Solve,  only: Equil_Solve
  use M_Equil_Tools,  only: Equil_FasTrace_EnTete,Equil_FasTrace_Write
  !
  use M_Global_Vars,  only: vSpc,vFas
  use M_Global_Vars,  only: vMixModel
  use M_Basis_Vars,   only: nCi,nAs,isW
  use M_Basis_Vars,   only: tNuFas,vOrdPr
  use M_Equil_Vars,   only: cEquMode,dXi
  use M_Equil_Vars,   only: vYesList
  use M_Equil_Vars,   only: vDeltaG_Fs,vAffScale,vLnAct,vFasMole
  use M_Equil_Vars,   only: vEquFas,nEquFas,T_EquPhase
  use M_Equil_Vars,   only: iCountEq,fTrcEq
  !
  logical,intent(in) :: NoMixture
  integer,intent(out):: iErr
  !
  integer :: nCp,nFs,nPur,nSol
  integer :: iFs,iMix,iPur
  integer :: FF
  integer :: nCheck,I,J,K
  integer :: nMix,nP
  integer :: iPrecip,iElimin
  integer :: nIts,iDo0,iDo1
  ! real(dp):: rDum= One
  ! integer :: vIPole(MaxPole)
  logical :: dXiTooLow
  logical :: Use_nIts= .true.  ! .false. !

  logical, allocatable:: vYesFas(:) ! phase is in current equil.phase list
  real(dp),allocatable:: vFasAff0(:)

  type(T_EquPhase),allocatable:: vEquFas0(:)
  type(T_MixPhase),allocatable:: vMixFas0(:)
  type(T_Phase),   allocatable:: vFas0(:)
  type(T_EquPhase):: EquFasNew
  ! type(T_MixPhase):: MixFas

  !~ type:: T_EquPhase
    !~ character(len=23):: NamEq
    !~ integer :: iSpc
    !~ integer :: iMix
    !~ real(dp):: vXPole(1:MaxPole)
    !~ !real(dp):: Grt
    !~ real(dp):: Mole
  !~ end type T_EquPhase

  if(iDebug>0) write(fTrc,'(/,A)') "< Equil_Eq1_Mix"

  FF=0
  if(iDebug>1 .and. (FF==0)) then
    !~ iCountEq=0
    !call Equil_FasTrace_EnTete(vFas,vYesList,FF)
    call GetUnit(FF)
    open(FF,file=trim(DirOut)//"_equil.log")
  end if

  cEquMode="EQ1" !-> flag in Equil_Solve
  dXi=       One !1.D3 !1.D-6 !reaction increment
  dXiTooLow= .false.

  nCheck= 0

  nCp= size(vOrdPr)
  !! nPur= count(vFas(:)%Typ=="PURE")
  nPur= count(vFas(:)%iSpc/=0)
  nMix= size(vMixModel)
  nSol= 1 !! size(vSolFas)

  if(NoMixture) nMix= 0

  nFs= nPur + nMix + nSol

  allocate(vYesFas(nFs)) ; vYesFas(:)= .false.
  allocate(vFasAff0(nFs))
  allocate(vEquFas0(nCp))

  if(nMix>0) then

    allocate(vFas0(nPur))
    vFas0(:)= vFas(1:nPur)
    deallocate(vFas)
    allocate(vFas(nFs))
    vFas(1:nPur)= vFas0(1:nPur)
    do K=1,nMix
      vFas(nPur+K)%NamFs= vMixModel(K)%Name
    end do
    deallocate(vFas0)

    allocate(vMixFas0(nMix))
    do K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    end do
    !
    allocate(tIPole(nMix,MaxPole))
    call MixModel_Init(nPur,vFas,vMixModel)

    !~ call MixPhase_XPole_Init(vMixModel,vMixFas_Xpole_Init)
    !
    !~ allocate(vDeltaGPole(nPur))
    !~ do iFs= 1,nPur
      !~ vDeltaGPole(iFs)= vFas(iFs)%Grt &
      !~ &               - dot_product(tNuFas(iFs,:),vSpc(vOrdPr(:))%G0rt)
    !~ end do

  end if
  !
  I=0
  do iPur=1,nPur
    if(vFas(iPur)%MolFs>FasMinim) then
      vYesFas(iPur)= .true.
      I=I+1
      vEquFas(I)%iPur= iPur
      vEquFas(I)%NamEq= trim(vFas(iPur)%NamFs)
      vEquFas(I)%iMix= 0
      vEquFas(I)%MolFs= vFas(iPur)%MolFs
    end if
  end do
  nEquFas= I
  !
  !-------------------------------------- loop on mineral saturation, --
  !------------------------ exit when all minerals have vPurYes false --
  iDo0= 0
  Do0: do
    iDo0= iDo0 +1

    iDo1= 0
    !---------------------- solve the system until no negative phases --
    Do1: do
      iDo1= iDo1 +1

      !~ if(iDebug>2) then
        !~ write(6,'(A)') "<-- CURRENT --"
        !~ do I=1,nEquFas
          !~ write(6,'(2X,A)') vEquFas(I)%NamEq
        !~ end do
        !~ write(6,'(A)') "</- CURRENT --"
      !~ end if

      if(nEquFas>0) then
        call EquFas_Clean
        call EquFas_Alloc(nCp,nEquFas)
        call EquFas_Update(nCp,nEquFas,vEquFas)
      end if

      !-------------------------------------------------- loop on dXi --
      Do2: do

        !------------------------------------------- solve the system --
        call Equil_Solve(nIts,iErr)
        !------------------------------------------/ solve the system --

        if(iErr==0) exit Do2 !when system solved-> exit

        !if necessary, try smaller dXi
        dXi= dXi /2.0_dp  !else do again with smaller step
        if(dXi<dXi_Minim) dXiTooLow=.true.
        if(dXiTooLow) exit Do0 !-> zannennagara,  found no solution ...

      end do Do2
      !--------------------------------------------------/loop on dXi --
      !
      if(nEquFas==0) exit Do1
      !
      if(Use_Nits) then
        if(nIts>12 .and. dXi<dXi_Minim) dXi= dXi /2.0_dp
        if(nIts<6  .and. dXi<dXi_Maxim) dXi= dXi *2.0_dp
      end if
      !
      if(iDebug>3) call Trace_1
      !
      !----------------------------------- exclude "near zero" phases --
      vEquFas0= vEquFas
      J= 0
      do I=1,nEquFas
        if(ABS(vEquFas0(I)%MolFs)>FasMinim) then
          J=J+1
          vEquFas(J)= vEquFas0(I)
          if(vEquFas0(I)%iPur>0) vYesFas(vEquFas0(I)%iPur)= .true.
          if(vEquFas0(I)%iMix>0) vYesFas(nPur+vEquFas0(I)%iPur)= .true.
        end if
      end do
      nEquFas= J

      vYesFas(:)= .false.
      vYesFas(vEquFas(1:nEquFas)%iPur)= .true.
      !----------------------------------/ exclude "near zero" phases --
      !
      !------------------------ check for phases with negative amount --
      if(nEquFas>0) then
        !
        !------------------------------------------------------ trace --
        ! if(iDebug>2) then
        !   write(6,'(A)') "<-- CURRENT --"
        !   do I=1,nEquFas
        !     write(6,'(G15.6,1X,A15,1X)',advance="NO") &
        !     & vEquFas(I)%Mole,trim(vEquFas(I)%NamEq)
        !     if(vEquFas(I)%iMix>0) then
        !       do J=1,vEquFas(I)%nPole
        !         write(6,'(G15.6,1X)',advance="NO") vEquFas(I)%vXPole(J)
        !       end do
        !     end if
        !     write(6,*)
        !   end do
        !   write(6,'(A)') "</- CURRENT --"
        ! end if
        !------------------------------------------------------/trace --
        !
        if(MINVAL(vEquFas(1:nEquFas)%MolFs)<-FasMinim) then
          iElimin= iMinLoc_R(vEquFas(1:nEquFas)%MolFs)
          if(iDebug>2) then
            do I=1,nEquFas
              write(6,'(2X,A)') vEquFas(I)%NamEq
            end do
            write(6,'(16X,2A)') "NEG= ", trim(vEquFas(iElimin)%NamEq)
          end if
          call EquFas_Remove(iElimin)

          if(vEquFas(iElimin)%iPur>0) &
          & vYesFas(vEquFas(iElimin)%iPur)= .false.
          if(vEquFas(iElimin)%iMix>0) &
          & vYesFas(nPur+vEquFas(iElimin)%iMix)= .false.

          cycle Do1 !--======< remove most "negative" phase and cycle ==
        end if
        !
      end if
      !
      exit Do1
      !-----------------------/ check for phases with negative amount --
      !
    end do Do1
    !---------------------/ solve the system until no negative phases --
    !
    if(iErr<0) exit Do0 != error in Newton -> exit
    !
    vFasAff0= Zero
    do iPur=1,nPur
      if(vYesList(iPur)) &
      & vFasAff0(iPur)= &
      & (  vDeltaG_Fs(iPur) &
      &  - dot_product(tNuFas(iPur,1:nCp),vLnAct(vOrdPr(1:nCp))) ) &
      & /vAffScale(iPur) !scaling ala EQ3/6
    end do
    !
    if(nMix>0) then
      !!call Mixture_Minimize(vFasAff0,vMixModel,vMixFas0)
      do iMix=1,nMix
        !vIPole(:)= tIPole(iMix,:)
        call Mixture_Minimize( &
        & vFasAff0,tIPole(iMix,:),vMixModel(iMix), &
        & vMixFas0(iMix))
      end do
      vFasAff0(nPur+1:nPur+nMix)= vMixFas0(1:nMix)%Grt
    end if
    !
    if(iDebug>3) call Trace_2
    !
    !---------------------------------------------------------- trace --
    !~ if(iDebug>2) then
      !~ write(6,'(A)') "<-- OVERSAT --"
      !~ do I=1,size(vFasAff0)
        !~ if (vFasAff0(I)<AffinIota) then
          !~ if(I<nPur) then
            !~ write(6,'(G15.6,1X,A15,1X)') vFasAff0(I),vFas(I)%NamFs
          !~ else
            !~ write(6,'(G15.6,1X,A15,1X)') vFasAff0(I),vFas(I)%NamFs
          !~ end if
        !~ end if
      !~ end do
      !~ write(6,'(A)') "</- OVERSAT --"
    !~ end if
    !----------------------------------------------------------/trace --
    !
    !----------------- add phase with highest supersaturation, if any --
    if(MINVAL(vFasAff0)<-AffinIota) then
      iPrecip= iMinLoc_R(vFasAff0)
      !-> the saturated phase with lowest vFasAff
      if(vYesFas(iPrecip)) iPrecip= 0
    else
      iPrecip= 0 != there is no saturated phase ==
      exit Do0 ! no phase found with affinity below -AffinIota > exit ==
    end if

    if(iPrecip>0) then

      if(iPrecip<=nPur) then

        vYesFas(iPrecip)= .true.

        EquFasNew%iPur= iPrecip
        EquFasNew%NamEq= trim(vFas(iPrecip)%NamFs)
        EquFasNew%iMix= 0
        EquFasNew%MolFs= Zero

      else

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

      !~ call CheckAssemblage(FF,nEquFas,EquFasNew,vEquFas,iElimin)
      !
      !~ if(iElimin>0) then
        !~ iPur= vEquFas(iElimin)%iPur
        !~ iMix= vEquFas(iElimin)%iMix
        !~ if(iPur > 0) vFas(iPur)%Mole= Zero
        !~ if(iMix > 0) vFas(nPur+iMix)%Mole= Zero
        !~ if(iDebug>2) &
        !~ & write(6,'(16X,2A)') "OUT= ",trim(vEquFas(iElimin)%NamEq)
        !~ call EquFas_Remove(iElimin)
      !~ end if

      if(iDebug>2) then
        do I=1,nEquFas
          write(6,'(2X,A)') vEquFas(I)%NamEq
        end do
        write(6,'(16X,2A)') "ADD= ",trim(EquFasNew%NamEq)
      end if

      if(nEquFas < size(vEquFas)) then
        vEquFas(nEquFas+1)= EquFasNew
        nEquFas= nEquFas +1
      end if

    end if

    if(iDebug>3) then !--=====================================< trace ==
      do I=1,nEquFas
        iFs= vEquFas(I)%iPur
        if(iFs>0) then
          write(6,'(A,I4,1X,A15,2(A,G15.6))') &
          & "Do0",iDo0, &
          & vFas(iFs)%NamFs, &
          & "/ lQsK=", -vFasAff0(iFs)/Ln10 ,&
          & "/ Mole=",vEquFas(I)%MolFs
        end if
      end do
      if(iDebug==4) call Pause_
    end if !==================================================</ trace ==

    !----------------/ add phase with highest supersaturation, if any --

    if(nEquFas==0) exit Do0

    !~ if(iDebug>2) pause
  end do Do0
  !-------------------------------------/ loop on mineral saturation, --
  !
  !--- save mineral mole numbers in vFas
  !~ call Save_EquPhase
  !
  deallocate(vFasMole)
  allocate(vFasMole(nPur+nMix))
  vFasMole(:)= Zero
  do I=1,nEquFas
    if(vEquFas(I)%iPur>0) vFasMole(vEquFas(I)%iPur)= vEquFas(I)%MolFs
    if(vEquFas(I)%iMix>0) vFasMole(vEquFas(I)%iMix+nPur)= vEquFas(I)%MolFs
  end do
  !~ pause
  !---/
  !
  call EquFas_Clean
  !
  !~ deallocate(vNewFas)
  deallocate(vFasAff0)
  deallocate(vEquFas0)
  if(nMix>0) then
    deallocate(vMixFas0)
    deallocate(tIPole)
  end if
  !
  if(FF>0) close(FF)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Equil_Eq1_Mix"
  !
contains

subroutine EquFas_Remove(N)
  integer,intent(in):: N
  integer:: I,J
  vEquFas0= vEquFas
  J= 0
  do I=1,nEquFas
    if(I/=N) then
      J=J+1
      vEquFas(J)= vEquFas0(I)
    end if
  end do
  nEquFas= J
end subroutine EquFas_Remove

subroutine Trace_1
  write(6,'(A)') "_"
  do I=1,nEquFas
    if(vEquFas(I)%iPur>0) &
    write(6,'(A,2(A,G15.6))') &
    & vFas(vEquFas(I)%iPur)%NamFs, &
    & "/ lQsK=",-vFasAff0(vEquFas(I)%iPur)/Ln10, &
    & "/ Mole=",vEquFas(I)%MolFs
  end do
end subroutine Trace_1

subroutine Trace_2
  write(6,'(6X,A)') "Do0"
  do I=1,nEquFas
    if(vEquFas(I)%iPur>0) &
    write(6,'(A,2(A,G15.6))') &
    & vFas(vEquFas(I)%iPur)%NamFs, &
    & "/ lQsK=",-vFasAff0(vEquFas(I)%iPur)/Ln10, &
    & "/ Mole=",vEquFas(I)%MolFs
  end do
end subroutine Trace_2

end subroutine Equil_Eq1

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

    ! find the indices of end-members of mixture vMixFas(iMix)
    ! in the pure phase list vFas(1:nPur)
    ! -> indices in tIPole(iMix,:)
    do P=1,nP
      P1= vMixModel(iModel)%vIPole(P)
      !-> index in vSpc
      !-> find the pure phases that points to this species
      P2= 0
      do iPur=1,nPur
        if(vFas(iPur)%iSpc==P1) then
          P2= iPur
          exit
        end if
      end do
      tIPole(iModel,P)= P2
      !
    end do

  end do

  return
end subroutine MixModel_Init

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
  use M_T_MixModel, only: MaxPole
  use M_Global_Vars,only: vMixModel
  use M_Basis_Vars, only: tAlfFs,tNuFas
  use M_Equil_Vars, only: vDeltaG_Eq,tAlfEq,tNuEq
  use M_Equil_Vars, only: vDeltaG_Fs,T_EquPhase
  !
  integer,intent(in):: nCp
  integer,intent(in):: nEquFas
  type(T_EquPhase),intent(in):: vEquFas(:)
  integer:: iPur,iMix,I,C
  integer:: nP
  integer:: vIPole(MaxPole)
  !
  do I=1,nEquFas
    !
    iPur= vEquFas(I)%iPur
    iMix= vEquFas(I)%iMix
    !
    if(iPur>0) then
      vDeltaG_Eq(I)= vDeltaG_Fs(iPur)
      tAlfEq(:,I)= tAlfFs(:,iPur)
      tNuEq(I,:)= tNuFas(iPur,:)
    end if
    !
    if(iMix>0) then
      nP= vEquFas(I)%NPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      vDeltaG_Eq(I)= SUM( vEquFas(I)%vXPole(1:nP) &
      &                 * vDeltaG_Fs(vIPole(1:nP)) ) &
      &            + SUM( vEquFas(I)%vXPole(1:nP) &
      &                 * vEquFas(I)%vLnAct(1:nP) )
      do C=1,nCp
        tAlfEq(C,I)= SUM( vEquFas(I)%vXPole(1:nP) &
        &               * tAlfFs(C,vIPole(1:nP)) )
        tNuEq(I,C)= SUM( vEquFas(I)%vXPole(1:nP) &
        &              * tNuFas(vIPole(1:nP),C) )
      end do
    end if
    !
  end do
  !
  write(sFMT,'(a,i3,a)') '(A,1X,',nCp,'(G12.3,1X))'
  do I=1,nEquFas
    write(73,'(A,G15.6,1X,A)') "DeltaG_Eq= ",vDeltaG_Eq(I),trim(vEquFas(I)%NamEq)
    write(73,sFMT)        "tAlfEq=    ",(tAlfEq(C,I),C=1,nCp)
    write(73,sFMT)        "tNuEq=     ",(tNuEq(I,C),C=1,nCp)
  end do
  write(73,'(A)') "===================================================="
  !
  return
end subroutine EquFas_Update

subroutine EquFasMix_Update(nCp,I,EquFas)
  use M_T_MixModel, only: MaxPole
  use M_Global_Vars,only: vMixModel
  use M_Basis_Vars, only: tAlfFs,tNuFas
  use M_Equil_Vars, only: vDeltaG_Eq,tAlfEq,tNuEq
  use M_Equil_Vars, only: vDeltaG_Fs,T_EquPhase
  !
  integer,intent(in):: nCp
  integer,intent(in):: I
  type(T_EquPhase),intent(in):: EquFas
  !
  integer:: C
  integer:: nP
  integer:: vIPole(MaxPole)
  !
  nP= EquFas%NPole
  vIPole(1:nP)= EquFas%vIPole(1:nP)
  vDeltaG_Eq(I)= SUM( EquFas%vXPole(1:nP) &
  &                 * vDeltaG_Fs(vIPole(1:nP)) ) &
  &            + SUM( EquFas%vXPole(1:nP) &
  &                 * EquFas%vLnAct(1:nP) )
  do C=1,nCp
    tAlfEq(C,I)= SUM( EquFas%vXPole(1:nP) &
    &               * tAlfFs(C,vIPole(1:nP)) )
    tNuEq(I,C)= SUM( EquFas%vXPole(1:nP) &
    &              * tNuFas(vIPole(1:nP),C) )
  end do
  !
  return
end subroutine EquFasMix_Update

subroutine EquFas_Clean
  use M_Equil_Vars,only: vDeltaG_Eq,tAlfEq,tNuEq
  
  if(allocated(vDeltaG_Eq)) deallocate(vDeltaG_Eq)
  if(allocated(tAlfEq))     deallocate(tAlfEq)
  if(allocated(tNuEq))      deallocate(tNuEq)
  
  return
end subroutine EquFas_Clean

subroutine Save_EquPhase
  use M_Equil_Vars, only: T_EquPhase
  use M_Global_Vars,only: vMixModel,vMixFas,vFas
  use M_Equil_Vars, only: nEquFas,vEquFas
  !
  !~ integer,intent(in):: nEquFas
  !~ type(T_EquPhase),intent(in):: vEquFas(:)
  !
  integer:: nPur,I,K,nP
  
  nPur= count(vFas(:)%iSpc>0)
  !
  vFas(:)%MolFs= Zero
  do I=1,nEquFas
    if(vEquFas(I)%iPur>0) vFas(vEquFas(I)%iPur)%MolFs= vEquFas(I)%MolFs
    if(vEquFas(I)%iMix>0) vFas(vEquFas(I)%iMix+nPur)%MolFs= vEquFas(I)%MolFs
  end do
  !
  K=0
  do I=1,nEquFas
    if(vEquFas(I)%iMix>0) K= K+1
  end do
  deallocate(vMixFas)
  allocate(vMixFas(K))
  !
  K=0
  if(size(vMixFas)>0) then
    do I=1,nEquFas
      if(vEquFas(I)%iMix>0) then
        K= K+1
        !
        nP= vMixModel(vEquFas(I)%iMix)%nPole
        vMixFas(K)%Name= trim(vEquFas(I)%NamEq)
        vMixFas(K)%iModel= vEquFas(I)%iMix
        vMixFas(K)%vXPole(1:nP)= vEquFas(I)%vXPole(1:nP)
        !
        vFas(nPur+K)%NamFs= trim(vEquFas(I)%NamEq)
        !! vFas(nPur+K)%Typ=   "MIXT"
        vFas(nPur+K)%iSpc=  0
        vFas(nPur+K)%iSol=  0
        vFas(nPur+K)%iMix=  K
        !
      end if
    end do
  end if
  
  return
end subroutine Save_EquPhase

subroutine Mixture_Minimize( &
& vAffPole,vIPole,MixModel, &
& MixFas) !,nFmix)
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
  use M_System_Vars,   only: TdgK,Pbar
  !
  use M_T_Phase,   only: T_Phase
  use M_T_MixPhase,only: T_MixPhase
  use M_T_MixModel,only: T_MixModel,MaxPole,Mix_Molecular
  !
  use M_Optimsolver_Theriak
  use M_MixModel_Optim
  !
  use M_Global_Vars,only: vFas
  !
  real(dp),        intent(in)   :: vAffPole(:)
  integer,         intent(in)   :: vIPole(:)
  type(T_MixModel),intent(in)   :: MixModel
  type(T_MixPhase),intent(inout):: MixFas
  !
  type(T_MixModel):: MM
  real(dp):: vX(MaxPole),vMu(MaxPole),vLnAct(MaxPole)
  real(dp):: vXmin(MaxPole),vMuMin(MaxPole)
  real(dp):: G,Gmin
  integer :: nP,K,P
  integer :: FF
  integer :: Multi
  real(dp):: S
  !
  real(dp):: TolX,DeltaInit
  integer :: its,nCallG
  logical :: OkConverge
  !
  FF= 0 ! log file
  TolX= MixMinim_TolX
  DeltaInit= 0.05D0
  !
  !! do I=1,size(vMixFas0)
    !
    !! iModel= vMixFas0(I)%iModel
    !! iModel= MixFas%iModel
    !! MM= vMixModel(iModel)
    MM= MixModel
    nP= MM%NPole
    !! vIPole(1:nP)= tIPole(iModel,1:nP)
    Multi= MM%vMulti(1)
    !
    !~ !---- compute affin' of all reactions between end-memb' and water --
    !~ do P= 1,nP
      !~ ! index of pole P in pure phase list
      !~ ! K= tIPole(iModel,P)
      !~ K= MM%vIPole(P)
      !~ ! compute affinity of reaction (pole P <> water )
      !~ !vAffPole(P)= vFas(K)%Grt &
      !~ !&          - dot_product(tNuFas(K,:), vSpc(vOrdPr(:))%G0rt) &
      !~ !&          - dot_product(tNuFas(K,:), vLnAct(vOrdPr(:)))
      !~ vAffPole(P)= vDeltaGPole(K) &
      !~ &          - dot_product(tNuFas(K,:), vLnAct(vOrdPr(:)))
    !~ end do
    !~ !---
    !
    if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
      !--------------------------------------------- analytic minimum --
      !
      Multi=  MM%vMulti(1)
      !
      S= Zero
      do P=1,nP
        ! K= MM%vIPole(P)
        ! vXmin(P)= FSafe_Exp(-vFasPur(MM%vIPole(P))%Grt /real(Multi))
        vXmin(P)= FSafe_Exp(-vAffPole(vIPole(P)) /real(Multi))
        S= S + vXmin(P)
      end do
      vXmin(1:nP)=  vXmin(1:nP) /S
      vLnAct(1:nP)= Multi*log(vXmin(1:nP))
      vMuMin(1:nP)= vAffPole(vIPole(1:nP)) + vLnAct(1:nP)
      do P=1,nP
        write(74,'(A,2G15.6,1X,A)') &
        & "vAffPole(P),vXmin(P)= ",&
        & vAffPole(vIPole(P)), &
        & vXmin(P), &
        & trim(vFas(vIPole(P))%NamFs)
      end do
      !pause
      Gmin= SUM(vXmin(1:nP) * vMuMin(1:nP))
      !
    else
      !----------------------------------------- numerical minimum(s) --
      !
      call MixModel_Optim_SetParams(TdgK,Pbar,MM)
      allocate(Mixmodel_Optim_vMu0rt(nP))
      allocate(Mixmodel_Optim_vLPole(nP))
      ! Mixmodel_Optim_vMu0rt(1:nP)= vFasPur(MM%vIPole(1:nP))%Grt
      Mixmodel_Optim_vMu0rt(1:nP)= vAffPole(vIPole(1:nP))
      Mixmodel_Optim_vLPole(1:nP)= .true.
      !
      Gmin= 1.0D30
      !
      do P=1,nP
        !
        !--- start from compos'nP close to end-member P
        vX(P)= One - 1.0D-3
        do K=1,nP
          if(K/=P) vX(K)= 1.0D-3/real(nP-1)
        end do
        !
        !~ vX(1:nP)= vMixFas_Xpole_Init(I)%tXPole(P,1:nP)
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
        !vMixFas_Xpole_Init(I)%tXPole(P,1:nP)= vX(1:nP)
        !
        if(G<Gmin) then
          Gmin= G
          vXmin(1:nP)= vX(1:nP)
          vMuMin(1:nP)= vMu(1:nP)
        end if
        !
      end do
      !pause
      !
      deallocate(Mixmodel_Optim_vMu0rt)
      deallocate(Mixmodel_Optim_vLPole)
      !
    end if
    !
    !!vMixFas0(i)%vLPole(:)= .true.
    !!vMixFas0(i)%vXPole(:)= vXmin(:)
    !!vMixFas0(i)%Grt= Gmin
    !!vMixFas0(i)%vLnAct(1:nP)= vMuMin(1:nP) -vAffPole(vIPole(1:nP))
    !
    MixFas%vLPole(:)= .true.
    MixFas%vXPole(:)= vXmin(:)
    MixFas%Grt= Gmin
    MixFas%vLnAct(1:nP)= vMuMin(1:nP) -vAffPole(vIPole(1:nP))
    !
  !! end do
  !
  !~ !-------------------------------------------------------- log files --
  !~ if(F1>0) then
    !~ write(F1,'(/,A,/)') "Mixture_Minimize,EndMember,Xi,Gi"
    !~ do I=1,size(vMixModel)
      !~ !if(vMixFas0(I)%Grt > 1.D-3) cycle
      !~ MM= vMixModel(I)
      !~ write(F1,'(2A)') "MODEL= ",MM%Name
      !~ do K=1,MM%NPole
        !~ write(F1,'(A,1X,G15.6)') &
        !~ & MM%vNamPole(K),        &
        !~ & vMixFas0(I)%vXPole(K)
      !~ end do
      !~ write(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
      !~ write(F1,*)
    !~ end do
    !~ write(F1,*)
  !~ end if
  !~ !
  !~ if(F2>0) then
    !~ do I=1,size(vMixModel)
    !~ !if(vMixFas0(I)%Grt < Zero) then
      !~ MM= vMixModel(I)
      !~ !
      !~ write(F2,'(A,A1)',advance="NO") MM%Name,T_
      !~ do K=1,MM%NPole
        !~ write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%vXPole(K),T_
      !~ end do
      !~ !do K=1,MM%NPole
      !~ !  write(F2,'(G15.6,A1)',advance="NO") vFasPur(MM%vIPole(K))%Grt,T_
      !~ !end do
      !~ write(F2,'(G15.6,A1)',advance="NO") vMixFas0(I)%Grt,T_
      !~ !
    !~ !end if
    !~ end do
    !~ write(F2,*)
  !~ end if
  !~ !-------------------------------------------------------/ log files --
  !
end subroutine Mixture_Minimize

end module M_Equil_1

