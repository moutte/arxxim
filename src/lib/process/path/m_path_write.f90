module M_Path_Write
  use M_Kinds
  use M_Trace,only: T_,fHtm,Pause_ !Stop_,fTrc,T_,iDebug
  
  implicit none
  
  private
  
  public:: Path_Write_Line
  public:: Path_Write_FasAff
  public:: Path_Write_FasEnTete
  public:: Path_Write_Distrib
  public:: Path_Files_Close
  
  integer::     & !
  & fQsK=    0, & !saturation states
  & fActiv=  0, & !log10(activities)
  & fPoten=  0, & !potentials /RT, Mu/RT(:)= Mu°(:)/RT + ln(act(:))
  & fMoles=  0, & !mole numbers
  & fMolal=  0, & !molality solutes
  & fGamma=  0    !gammas
  !
  logical:: bOpenDist=.false.
  !
  character(len=15),allocatable:: vNamFasYes(:)

contains

subroutine Path_Write_FasEnTete
!--
!-- write entete for the file of affinities of pure phases
!-- (i.e. non-aqu'species)
!--
  use M_IOTools,   only: GetUnit
  use M_Files,     only: cTitle,DirOut, Files_Index_Write
  use M_Global_Vars,only: vFas
  !
  integer :: iFs
  !
  call GetUnit(fQsK)
  open(fQsK,file=trim(DirOut)//"_minqsk.restab")
  call Files_Index_Write(fHtm, &
  & trim(DirOut)//"_minqsk.restab", &
  & "Q/K of all phases along path ")
  !
  write(fQsK,'(A,A1)',advance= 'NO') "Count",T_ !,"pH",T_
  do iFs=1,size(vFas)
    write(fQsK,'(A,A1)',advance= 'NO') trim(vFas(iFs)%NamFs),T_
  end do
  write(fQsK,*)
  !
end subroutine Path_Write_FasEnTete

subroutine Path_Write_FasAff(iStep)
!--
!-- write affinity of pure phases (i.e. non-aqu'species)
!--
  use M_Numeric_Const,only: Ln10
  use M_Global_Vars,  only: vFas,vSpc
  use M_Basis_Vars,   only: tNuFas,vOrdPr
  !
  integer,intent(in):: iStep
  !
  integer :: iFs
  real(dp):: X
  !
  write(fQsK,'(I4,A1)',advance= 'NO') iStep,T_
  do iFs=1,size(vFas)
    X= &
    & vFas(iFs)%Grt &
    - dot_product( tNuFas(iFs,:), &
    &              vSpc(vOrdPr(:))%Dat%LAct+vSpc(vOrdPr(:))%G0rt )
    X= - X /SUM(ABS(tNuFas(iFs,:))) /Ln10
    !!write(fQsK,'(G15.6,A1)',advance="no") -vFasAff(iFs) /Ln10,T_ !-> log10(QsK)= - Affinity /Ln10
    write(fQsK,'(G15.6,A1)',advance="no") X,T_ !-> log10(QsK)= - Affinity /Ln10
  end do
  write(fQsK,*)
  !
end subroutine Path_Write_FasAff

subroutine Path_Write_Distrib(N)
!--- distribution of all elements among their different species --
!--- to be used after Equil_Save (to update vSpc%Mole) --
!--
!-- Element iEl is distributed among all Species iAq that have
!-- tFormula(iEl,iAq)/=0
!--
!-- if Speciation succeeds, 
!-- then we should have
!--   vCpn(iEl)%Mole=dot_product(tFormula(iEl,1:nAq),vSpc(1:nAq)%Dat%Mole)),
!-- except for element iBal or for a mobile element
!--
!-- Output format is designed for easy data processing w/ spreadsheet
!-- (sort according to column vNamEl, -> produce tables for the different elements)
!-- currently no special treatment for Redox !!! should do something...
!--
  use M_IoTools,only: GetUnit
  use M_Numeric_Const,  only: Ln10
  use M_SolModel_Tools,only: Solmodel_pHpE
  !
  use M_Global_Vars, only: vEle,nAq,vSpc,tFormula
  use M_System_Vars, only: vCpn
  use M_Files,       only: DirOut
  use M_Basis_Vars,  only: isW,iOx,isH_,isO2
  !
  integer ::fDist,N,iAq,iEl,nCp
  real(dp)::Tot, X, pH_,pE_
  !
  call Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  !
  nCp=size(vCpn)
  !
  if(.not.bOpenDist) then
    call EnTeteDistrib
    bOpenDist=.true.
  end if
  !
  do iEl=2,nCp
    !
    call GetUnit(fDist)
    open(fDist, &
    & file=trim(DirOut)//"_distrib_"//trim(vEle(vCpn(iEl)%iEle)%NamEl)//".restab", &
    & ACCESS='APPend')
    !
    Tot=dot_product(tFormula(iEl,1:nAq),vSpc(1:nAq)%Dat%Mole)
    if(trim(vCpn(iEl)%Statut)/="INERT") vCpn(iEl)%Mole=Tot
    !
    write(fDist,'(I3,A1,A,A1,F7.3,A1,F7.3,A1)',advance="no") &
    &             N,T_,trim(vEle(vCpn(iEl)%iEle)%NamEl),T_,pH_,T_,pE_,T_
    !
    do iAq=2,nAq
      !if(iAq/=iH_) then
        if(tFormula(iEl,iAq)/=0) then
          X=vSpc(iAq)%Dat%Mole/tFormula(iEl,iAq)/Tot*1000.0D0
          !divides vMolF of species iAq by number of moles of iEl in one mole of iAq
          !if(X>1.0D-3) then; write(fDist,'(F12.3,A1)',advance="no") X,T_ 
          !             else; write(fDist,'(A12,A1)',  advance="no") "0",T_; end if 
          write(fDist,'(G15.6,A1)',advance="no") X,T_
        end if
      !end if
    end do
    !
    do iAq=2,nAq
      if(tFormula(iEl,iAq)/=0) &
      & write(fDist,'(G15.6,A1)',advance="no") 1.0E6*vSpc(iAq)%Dat%Mole,T_
      !*1.0E6 -> output in MICROMOLES !!!
    end do
    do iAq=2,nAq
      if(tFormula(iEl,iAq)/=0) &
      & write(fDist,'(G15.6,A1)',advance="no") exp(vSpc(iAq)%Dat%LGam),T_ !vSpc(iAq)%LnGam/Ln10
    end do
    !
    write(fDist,'(2(G15.6,A1))') 1.0E6*Tot,T_,1.0E6*vCpn(iEl)%Mole
    close(fDist)
    !
  end do
  !
end subroutine Path_Write_Distrib

subroutine EnteteDistrib
!--
!-- write head line for distribution of elements among species --
!--
  use M_IoTools,     only: GetUnit
  use M_Global_Vars, only: nAq,vEle,vSpc,tFormula
  use M_System_Vars, only: vCpn
  use M_Files,       only: cTitle,DirOut
  !
  integer::fDist,iAq,iEl
  !
  do iEl=2,size(vCpn)
    !
    call GetUnit(fDist)
    open(fDist,file=trim(DirOut)//"_distrib_"//trim(vEle(vCpn(iEl)%iEle)%NamEl//".restab"))
    !
    !if(iEl/=iH_) then
    write(fDist,'(A,/A)') &
    & "!."//trim(cTitle), &
    & "!.species distribution for a given element, relative (permil) and absolute (micromoles/kgH2O)"
    !
    write(fDist,'(4(A,A1))',advance="no") &
    &            "count",T_,trim(vEle(vCpn(iEl)%iEle)%NamEl),T_,"pH",T_,"pE",T_
    !
    do iAq=2,nAq
      if(tFormula(iEl,iAq)/=0) &
      &  write(fDist,'(A,A1)',  advance="no") trim(vSpc(iAq)%NamSp)//"_rel",T_
    end do
    do iAq=2,nAq
      if(tFormula(iEl,iAq)/=0) &
      &  write(fDist,'(A,A1)',  advance="no") trim(vSpc(iAq)%NamSp),T_
    end do
    do iAq=2,nAq
      if(tFormula(iEl,iAq)/=0) &
      &  write(fDist,'(A,A1)',  advance="no") trim(vSpc(iAq)%NamSp)//"_Gam",T_
    end do
    !
    write(fDist,'(A,A1,A)',advance="no") "TOT,umoles",T_,"Check"
    write(fDist,*)
    !
    !end if
    close(fDist)
  end do
  !
  return
end subroutine EnTeteDistrib

subroutine Path_Write_Line(WrCount,WrCod,vYes,WrTitl)
!--
!-- one result on one line
!-- writes on several files:
!-- species log(activity) on _activ,
!-- species molalities on _molal,
!-- element abundance on _moles
!--
  use M_IOTools
  use M_Files,        only: DirOut,cTitle,Files_Index_Write
  use M_Dtb_Const,    only: T_CK
  use M_Dtb_Const,    only: R_JK,Tref
  use M_Numeric_Const,only: Ln10
  use M_SolModel_Tools, only: Solmodel_pHpE,Solmodel_CalcMolal
  use M_T_Species,    only: Species_EntropyZero
  !
  use M_Global_Vars,only: vEle,vSpc,vFas,nAq,SolModel !!,Solvent
  use M_System_Vars,only: TdgK,Pbar,vCpn
  use M_Basis_Vars, only: isW,MWsv,iOx,isH_,isO2,tStoikio,tAlfFs
  use M_Basis_Vars, only: vLCi,vLAs,vLAx,vLMx !,nAx,nMx,vYesList
  use M_Path_Vars,  only: DimPath,tPathResults
  !
  integer,              intent(in):: WrCount
  character(len=3),     intent(in):: WrCod
  logical,     optional,intent(in):: vYes(:)
  character(*),optional,intent(in):: WrTitl
  !
  logical,parameter:: OutputInGrams = .false.  !! .true.
  !
  real(dp):: vMolal(1:nAq)
  real(dp):: ZBal,ZRdx,IonStr,xTotF
  real(dp):: pH_,pE_,X,nOx
  integer :: iPr,I,N,K,nCp
  !
  real(dp),allocatable:: vS0Ele(:)
  real(dp),allocatable:: vPot(:)
  
  allocate(vS0Ele(size(vSpc)))
  allocate(vPot(size(vSpc)))
  !
  nCp= size(vCpn)
  !
  do I=1,size(vSpc)
    !vS0Ele(I)= dot_product(vSpc(I)%vStoikio(1:N),vEle(1:N)%S0) /real(vSpc(I)%vStoikio(0))
    vS0Ele(I)= Species_EntropyZero(vEle,vSpc(I))
  end do
  !
  call Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  !
  !--------------------------------------------------ACTIV'COEFFS.HEADER
  if(fGamma==0) then
    !
    call GetUnit(fGamma)
    open(fGamma,file=trim(DirOut)//"_gamma.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_gamma.restab", &
    & "PATH: log10(activity coeff's)")
    !
    write(fGamma,'(A,A1)',   advance="no") "Count",T_
    write(fGamma,'(3(A,A1))',advance="no") "TempC",T_,"Pbar",T_,"RhoW",T_
    write(fGamma,'(3(A,A1))',advance="no") "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    if(iOx/=0) &
    & write(fGamma,'(A2,A1,A5,A1)',advance="no") "pE",T_,"RedOx",T_
    !
    !write component names
    do iPr=1,size(vCpn) 
      write(fGamma,'(A,A1)',advance="no") trim(vEle(vCpn(iPr)%iEle)%NamEl),T_
    end do
    !
    !write species names (aqu'species only)
    do I=1,size(vSpc)
      if(vSpc(I)%Typ(1:3)=="AQU") & !!(vLCi(I).or.vLAs(I).or.vLAx(I))
      & write(fGamma,'(A,A1)',advance="no") trim(vSpc(I)%NamSp),T_
    end do
    !
    write(fGamma,*)
    !
  end if
  !-------------------------------------------------/ACTIV'COEFFS.HEADER
  !
  !----------------------------------------------------MOLALITIES.HEADER
  if(fMolal==0) then
    !
    call GetUnit(fMolal)
    open(fMolal,file=trim(DirOut)//"_molal.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_molal.restab", &
    & "PATH: molalities")
    !
    write(fMolal,'(A,A1)',   advance="no") "Count",T_
    write(fMolal,'(3(A,A1))',advance="no") "TempC",T_,"Pbar",T_,"RhoW",T_
    write(fMolal,'(3(A,A1))',advance="no") "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    if(iOx/=0) &
    & write(fMolal,'(A2,A1,A5,A1)',advance="no") "pE",T_,"RedOx",T_
    !
    !--- write component names
    do iPr=1,size(vCpn) 
      write(fMolal,'(A,A1)',advance="no") trim(vEle(vCpn(iPr)%iEle)%NamEl),T_
    end do
    !
    !--- write species names (aqu'species only)
    do I=1,size(vSpc)
      if(vSpc(I)%Typ(1:3)=="AQU") & !!(vLCi(I).or.vLAs(I).or.vLAx(I))
      & write(fMolal,'(A,A1)',advance="no") trim(vSpc(I)%NamSp),T_
    end do
    !
    write(fMolal,*)
    !
  end if
  !---------------------------------------------------/MOLALITIES.HEADER
  !
  !----------------------------------------------------POTENTIALS.HEADER
  if(fPoten==0) then
    !
    call GetUnit(fPoten)
    open(fPoten,file=trim(DirOut)//"_potent.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_potent.restab", &
    & "PATH result: log10(species activities)")
    !
    if(present(WrTitl)) write(fPoten,'(A,A1)',advance="no") ".Title",T_
    write(fPoten,'(A,A1)',   advance="no") "Count",T_
    write(fPoten,'(3(A,A1))',advance="no") "TempC",T_,"Pbar",T_,"RhoW",T_
    write(fPoten,'(3(A,A1))',advance="no") "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    if(iOx/=0) &
    & write(fPoten,'(A2,A1,A5,A1)',advance="no") "pE",T_,"RedOx",T_
    !
    !--- write component names
    do iPr=1,size(vCpn) 
      write(fPoten,'(A,A1)',advance="no") trim(vEle(vCpn(iPr)%iEle)%NamEl),T_
    end do
    !
    !--- write species names (all aqu' + mobile non'aqu')
    do I=1,size(vSpc)
      if(vLCi(I).or.vLAs(I).or.vLAx(I).or.vLMx(I)) &
      & write(fPoten,'(A,A1)',advance="no") trim(vSpc(I)%NamSp),T_
    end do
    write(fPoten,*)
    !
    !-------------------------------------------------------------------
    !! !
    !! if(present(WrTitl)) write(fPoten,'(A,A1)',advance="no") ".Title",T_
    !! write(fPoten,'(A,A1)',   advance="no") "Count",T_
    !! write(fPoten,'(3(A,A1))',advance="no") "TempC",T_,"Pbar",T_,"RhoW",T_
    !! write(fPoten,'(3(A,A1))',advance="no") "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !! !
    !! if(iOx/=0) &
    !! & write(fPoten,'(A2,A1,A5,A1)',advance="no") "pE",T_,"RedOx",T_
    !! !
    !! !--- write component names
    !! do iPr=1,size(vCpn) 
    !!   write(fPoten,'(A,A1)',advance="no") trim(vEle(vCpn(iPr)%iEle)%NamEl),T_
    !! end do
    !! !
    !! do I=1,size(vSpc)
    !!   if(vLCi(I).or.vLAs(I).or.vLAx(I).or.vLMx(I)) &
    !!   & write(fPoten,'(G15.6,A1)',advance="no") vSpc(I)%G0rt,T_
    !! end do
    !! write(fPoten,*)
    !
  end if
  !--------------------------------------------------/ POTENTIALS.HEADER
  !
  !--------------------------------------------------- ACTIVITIES.HEADER
  if(fActiv==0) then
    !
    call GetUnit(fActiv)
    open(fActiv,file=trim(DirOut)//"_activ.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_activ.restab", &
    & "PATH result: log10(species activities)")
    !
    if(present(WrTitl)) write(fActiv,'(A,A1)',advance="no") ".Title",T_
    write(fActiv,'(A,A1)',   advance="no") "Count",T_
    write(fActiv,'(3(A,A1))',advance="no") "TempC",T_,"Pbar",T_,"RhoW",T_
    write(fActiv,'(3(A,A1))',advance="no") "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    if(iOx/=0) &
    & write(fActiv,'(A2,A1,A5,A1)',advance="no") "pE",T_,"RedOx",T_
    !
    !--- write component names
    do iPr=1,size(vCpn) 
      write(fActiv,'(A,A1)',advance="no") trim(vEle(vCpn(iPr)%iEle)%NamEl),T_
    end do
    !
    !--- write species names (all aqu' + mobile non'aqu')
    do I=1,size(vSpc)
      if(vLCi(I).or.vLAs(I).or.vLAx(I).or.vLMx(I)) &
      & write(fActiv,'(A,A1)',advance="no") trim(vSpc(I)%NamSp),T_
    end do
    write(fActiv,*)
    !
  end if
  !-------------------------------------------------/ACTIVITIES.HEADER--
  !
  !call OutStrVec(fMoles,vCpn(1:nCp)%Mole,I=WrCount,S="TotEle=",C="F")
  !_____________________________________Output
  !call OutStrVec(fMoles,vSpc(1:nAq)%Dat%Mole,I=WrCount,S="MolNum=",C="F")
  !
  ZBal=dot_product(vSpc(1:nAq)%Dat%Mole,vSpc(1:nAq)%Z)
  !
  call Solmodel_CalcMolal(vSpc,isW,vMolal,IonStr)
  !
  !--------------------------------------------------left part of tables
  !
  !-------------------------------------------------ACTIV'COEFF'.resultS
  write(fGamma,'(I3,A1)',      advance="no") WrCount,T_
  write(fGamma,'(3(G15.3,A1))',advance="no") TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  write(fGamma,'(3(G15.6,A1))',advance="no") IonStr,T_,pH_,T_,ZBal,T_
  !
  !---------------------------------------------------ACTIVITIES.resultS
  if(present(WrTitl)) write(fActiv,'(A,A1)',advance="no") trim(WrTitl),T_
  write(fActiv,'(I3,A1)',      advance="no") WrCount,T_
  write(fActiv,'(3(G15.3,A1))',advance="no") TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  write(fActiv,'(3(G15.6,A1))',advance="no") IonStr,T_,pH_,T_,ZBal,T_
  !
  !---------------------------------------------------POTENTIALS.resultS
  if(present(WrTitl)) write(fActiv,'(A,A1)',advance="no") trim(WrTitl),T_
  write(fPoten,'(I3,A1)',      advance="no") WrCount,T_
  write(fPoten,'(3(G15.3,A1))',advance="no") TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  write(fPoten,'(3(G15.6,A1))',advance="no") IonStr,T_,pH_,T_,ZBal,T_
  !
  !---------------------------------------------------MOLALITIES.resultS
  write(fMolal,'(I3,A1)',      advance="no") WrCount,T_
  write(fMolal,'(3(G15.3,A1))',advance="no") TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  write(fMolal,'(3(G15.6,A1))',advance="no") IonStr,T_,pH_,T_,ZBal,T_
  !
  if(iOx/=0) then
    ZRdx= dot_product(vSpc(1:nAq)%Dat%Mole,tStoikio(iOx,1:nAq))
    write(fActiv,'(2(G15.6,A1))',advance="no") pE_,T_,ZRdx,T_
    write(fPoten,'(2(G15.6,A1))',advance="no") pE_,T_,ZRdx,T_
    write(fGamma,'(2(G15.6,A1))',advance="no") pE_,T_,ZRdx,T_
    write(fMolal,'(2(G15.6,A1))',advance="no") pE_,T_,ZRdx,T_
  end if
  !
  do iPr=1,size(vCpn)
    !
    xTotF= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole)
    !
    write(fActiv,'(G15.8,A1)',advance="no") xTotF,T_
    write(fPoten,'(G15.8,A1)',advance="no") xTotF,T_
    write(fGamma,'(G15.8,A1)',advance="no") xTotF,T_
    write(fMolal,'(G15.8,A1)',advance="no") xTotF,T_
    !
  end do
  !-----------------------------------------------/left part of tables--
  !
  !!RdxBal=Zero
  !!if(iOx/=0) then !redox balance
  !!  RdxBal=                    dot_product(tAlfPr(iOx,1:nCi),        vMolF(1:        nCi))
  !!  if(nAs>0) RdxBal= RdxBal + dot_product(tAlfAs(iOx,1:    nAs),    vMolF(nCi+1:    nCi+nAs))
  !!  if(nAx>0) RdxBal= RdxBal + dot_product(tAlfPr(iOx,nCi+1:nCi+nAx),vMolF(nCi+nAs+1:nCi+nAs+nAx))
  !!end if
  !
  !
  !!! do I=1,size(vSpc)
  !!!   if(vLCi(I).or.vLAs(I).or.vLAx(I).or.vLMx(I)) &
  !!!   & write(fActiv,'(G15.6,A1)',advance="no") vSpc(I)%LnAct/Ln10,T_
  !!! end do
  !!! write(fActiv,*)
  !
  !-----------------------------------------------RIGHT PART OF TABLES--
  !-------------------------------------------------ACTIVITIES.resultS--
  !! print *,"vLMx=",size(vLCi),size(vLAs),size(vLAx),size(vLMx)
  !! do i=1,size(vLMx)
  !! print *,vLMx(i)
  !! end do
  !! print *,"OutStrVec"  ;  call pause_
  call OutStrVec( &
  & fActiv,vSpc(:)%Dat%LAct/Ln10) !, &
  ! & vLc(:)= vLCi &
  ! & .or.vLAs &
  ! & .or.vLAx &
  ! & .or.vLMx)
  !
  !---------------------------------------------------POTENTIALS.resultS
  vPot(:)= vSpc(:)%Dat%LAct +vSpc(:)%G0rt !-> "reduced"potential, dimensionless
  vPot(:)= vPot(:) *R_JK *TdgK            !-> potential in Joule
  ! vPot(:)= vPot(:) - Tref*vS0Ele(:)       !-> in Berman-Brown convention
  !
  call OutStrVec( &
  & fOut= fPoten,&
  & Vec=  vPot(:), &
  & vL=   vLCi.or.vLAs.or.vLAx.or.vLMx)
  !
  !--------------------------------------------------ACTIV'COEFF'.resultS
  call OutStrVec( &
  & fOut= fGamma,&
  & Vec=  vSpc(:)%Dat%LGam/Ln10, &
  & vL=   vLCi.or.vLAs.or.vLAx)
  !
  !---------------------------------------------------MOLALITIES.resultS
  call OutStrVec( &
  & fOut= fMolal,&
  & Vec=  vMolal(:))
  !------------------------------------------------/RIGHT PART OF TABLES
  !
  !---------------------------------------------------------------------
  !--------------------------------------------------MOLE.NUMBERS.HEADER
  if(fMoles==0) then
    !
    !if(allocated(tPathResults)) deallocate(tPathResults)
    !if(allocated(vNamFasYes)) deallocate(vNamFasYes)
    !
    call GetUnit(fMoles)
    open(fMoles,file=trim(DirOut)//"_moles.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_moles.restab", &
    & "PATH: mole nrs of elements, influid, in system, mole nrs of other phases")
    !
    if(present(WrTitl)) write(fMoles,'(A,A1)',advance="no") ".Title",T_
    write(fMoles,'(4(A,A1))',advance="no") "count",T_,"TdgC",T_,"Pbar",T_,"RhoW",T_
    write(fMoles,'(2(A,A1))',advance="no") "pH",T_,"pE",T_
    !
    !----------------------------------------------write component names
    do iPr=1,size(vCpn)
      write(fMoles,'(A,A1)',advance="no") trim(vEle(vCpn(iPr)%iEle)%NamEl),T_
    end do
    !
    !-----------------------------------------write non'aqu'phases names
    if(present(vYes) .and. WrCod(1:2)=="EQ") then
      !
      !--------------------------------------------allocate BUFFER TABLE
      !---------will be used to build a result file without zero columns
      !--------used especially for EQUPATH with large number of minerals
      !
      if(count(vYes)>0) then
        allocate(vNamFasYes(count(vYes)))
        allocate(tPathResults(DimPath,2+2*nCp+count(vYes)))
        tPathResults(:,:)= Zero
      end if
      !-------------------------------------------/allocate BUFFER TABLE
      !
      do iPr=1,size(vCpn)
        write(fMoles,'(A,A1)',advance="no") &
        & trim(vEle(vCpn(iPr)%iEle)%NamEl)//"Tot",T_
      end do
      !
      K=0
      do I=1,size(vFas)
        if(vYes(I)) then
          write(fMoles,'(A,A1)',advance="no") trim(vFas(I)%NamFs),T_
          K=K+1
          if(allocated(vNamFasYes)) vNamFasYes(K)= trim(vFas(I)%NamFs)
        end if
      end do
    end if
    !
    write(fMoles,*)
    !
  end if
  !-------------------------------------------------/MOLE.NUMBERS.HEADER
  !
  !-------------------------------------------------MOLE.NUMBERS.resultS
  !write(fMoles,'(A4,A1,I3,A1)',advance="no") "ELEM",T_,WrCount,T_
  if(present(WrTitl)) write(fMoles,'(A,A1)',advance="no") trim(WrTitl),T_
  !
  write(fMoles,'(I4,A1,3(G15.3,A1))',advance="no")  &
  & WrCount,T_,TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  !
  write(fMoles,'(2(G15.3,A1))',advance="no")  pH_,T_,pE_,T_
  !
  if(allocated(tPathResults)) then !---------------write IN BUFFER TABLE
    tPathResults(WrCount,1)= TdgK -T_CK
    tPathResults(WrCount,2)= Pbar
  end if
  !
  !-----------------------------------------write total element in fluid
  do iPr=1,size(vCpn)
    X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> mole nrs in fluid
    write(fMoles,'(G15.8,A1)',advance="no") X,T_
    if(allocated(tPathResults)) &
    & tPathResults(WrCount,iPr+2)= X !-------------write IN BUFFER TABLE
  end do
  !
  !---------------------------------write total element in fluid + min's
  if(present(vYes) .and. WrCod(1:2)=="EQ") then
    !
    do iPr=1,size(vCpn)
      !
      X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> mole nrs in fluid
      !
      do I=1,size(vFas)
        if(vYes(I) .and. vFas(I)%MolFs>Zero) &
        & X= X + tAlfFs(iPr,I) *vFas(I)%MolFs !-> mole nrs in fluid + min's
      end do
      !
      write(fMoles,'(G15.8,A1)',advance="no") X,T_
      !
      if(allocated(tPathResults)) &
      & tPathResults(WrCount,iPr+nCp+2)= X !-------write IN BUFFER TABLE
      !
    end do
    !
    K= nCp*2 +2
    do I=1,size(vFas) !mineral mole numbers from Equil_N
      !if(.not.vFasYes(I))  vFas(I)%MolFs=-999.D0 
      if(vYes(I)) then
        !
        K=K+1
        if(vFas(I)%MolFs>Zero) then
          !
          !! nOx= One !*tFormula(iFs,1)
          !! if(tAlfFs(1,I)/=0) nOx= tAlfFs(1,I) !-> nr of oxygens in formula
          !! write(fMoles,'(G15.6,A1)',advance="no") vFas(I)%MolFs*nOx,T_ !  &
          !
          !! write(fMoles,'(G15.6,A1)',advance="no") 
          !! & vFas(I)%MolFs /vCpn(1)%Mole/MWsv,T_
          !!-> mole nrs of phase / kgH2O
          !
          nOx= One !*tFormula(iFs,1)
          !
          if(OutputInGrams) then
            !
            write(fMoles,'(G15.6,A1)',advance="no") &
            & vFas(I)%MolFs*vFas(I)%WeitKg*1.0D3,T_
            !
            !-- write IN BUFFER TABLE --
            if(allocated(tPathResults)) &
            & tPathResults(WrCount,K)= vFas(I)%MolFs*vFas(I)%WeitKg*1.0D3
            !
          else
            !
            if(tAlfFs(1,I)/=0) nOx= tAlfFs(1,I) !-> nr of oxygens in formula
            write(fMoles,'(G15.6,A1)',advance="no") vFas(I)%MolFs*nOx,T_ !  &
            !
            !-- write IN BUFFER TABLE --
            if(allocated(tPathResults)) &
            & tPathResults(WrCount,K)= vFas(I)%MolFs*nOx
            !
          end if
          !
        else
          !
          write(fMoles,'(A,A1)',advance="no") "0.0",T_
          if(allocated(tPathResults)) &
          & tPathResults(WrCount,K)= Zero
          !
        end if
        !
        !
      end if
    end do
    !
  end if
  !--------------------------------/write total element in fluid + min's
  write(fMoles,*)
  !------------------------------------------------/MOLE.NUMBERS.resultS
  !---------------------------------------------------------------------
  !
  deallocate(vS0Ele)
  deallocate(vPot)
  !
  return
end subroutine Path_Write_Line

subroutine WriteBuffer
  use M_IOTools
  use M_Files,    only: DirOut
  !
  use M_Global_Vars,only: vFas,vEle
  use M_System_Vars,only: vCpn
  use M_Path_Vars,only: tPathResults
  !
  integer:: FF,I,J,K,nYes,nCp
  logical,allocatable:: vIsZero(:)
  !
  nCp= size(vCpn)
  nYes= size(tPathResults,2) -(2+2*nCp)
  allocate(vIsZero(nYes))
  !
  do J=1,nYes
    vIsZero(J)= (ABS(MAXVAL(tPathResults(:,2+2*nCp+J)))<1.D-16)
  end do
  !
  !print *,"=====nYes-count(vIsZero)===",nYes-count(vIsZero)
  !return
  !
  call GetUnit(FF)
  open(FF,file=trim(DirOut)//"_moles_clean.restab")
  !
  !-------------------------------------------------------------HEADER--
  write(FF,'(A,A1)',advance="no") "count",T_
  write(FF,'(2(A,A1))',advance="no") "TdgC",T_,"Pbar",T_
  !~ write(FF,'(3(A,A1))',advance="no") "RhoW",T_,"pH",T_,"pE",T_
  !~ !
  !-----------------------------------------------write component names--
  do I=1,size(vCpn)
    write(FF,'(A,A1)',advance="no") &
    & trim(vEle(vCpn(I)%iEle)%NamEl)//"Fluid",T_
  end do
  !
  do I=1,size(vCpn)
    write(FF,'(A,A1)',advance="no") &
    & trim(vEle(vCpn(I)%iEle)%NamEl)//"Tot",T_
  end do
  !-----------------------------------------write non'aqu'phases names--
  K=0
  do I=1,nYes
    K=K+1
    if(.not. vIsZero(I)) then
      write(FF,'(A,A1)',advance="no") trim(vNamFasYes(K)),T_
    end if
  end do
  write(FF,*)
  !------------------------------------------------------------/HEADER--
  !
  do I=1,size(tPathResults,1)
    !
    write(FF,'(I3,A1)',advance="no") I,T_
    write(FF,'(G15.8,A1)',advance="no") tPathResults(I,1),T_ !! TdgC
    write(FF,'(G15.8,A1)',advance="no") tPathResults(I,2),T_ !! Pbar
    !
    do J=1,2*nCp
      write(FF,'(G15.8,A1)',advance="no") tPathResults(I,J+2),T_
    end do
    !
    do J=1,nYes
      if(.not. vIsZero(J)) &
      & write(FF,'(G15.8,A1)',advance="no") tPathResults(I,J+2+2*nCp),T_
    end do
    !
    write(FF,*)
    !
  end do
   !
  deallocate(vIsZero)
  !
  close(FF)
  !
end subroutine WriteBuffer

subroutine Path_Files_Close
  use M_Path_Vars,only: tPathResults
  !
  if(fActiv>0) then  ;  write(fActiv,*)  ;  close(fActiv) ;  fActiv=0  ;  end if
  if(fPoten>0) then  ;  write(fPoten,*)  ;  close(fPoten) ;  fPoten=0  ;  end if
  if(fMolal>0) then  ;  write(fMolal,*)  ;  close(fMolal) ;  fMolal=0  ;  end if
  if(fGamma>0) then  ;                      close(fGamma) ;  fGamma=0  ;  end if
  if(fQsk>0)   then  ;                      close(fQsk)   ;  fQsk=0    ;  end if
  !
  if(fMoles>0) then
  
    write(fMoles,*)
    close(fMoles)
    fMoles=0
    
    if(allocated(tPathResults)) then
      call WriteBuffer
      deallocate(tPathResults)
      deallocate(vNamFasYes)
    end if
    
  end if
  !
end subroutine Path_Files_Close

end module M_Path_Write

