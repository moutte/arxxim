module M_Equil_Tools
!--
!-- common tools for equilibrium calculations
!--
  use M_Trace,only: iDebug,fTrc,fHtm,T_,Stop_,Pause_,Warning_
  use M_Info_Value
  use M_Kinds
  !
  implicit none
  !
  private
  !
  !!public:: Equil_Solve
  !
  public:: Equil_Zero      !for a new speciation or equilibrium
  public:: Equil_Restart
  public:: Equil_Save      !store results of Equil_Calc in vCpn, vSpc 
  public:: Equil_Errors
  public:: Equil_Errors_Warning
  !
  public:: Equil_CalcFasAff
  !
  public:: Equil_Trace_Init
  public:: Equil_Trace_Close
  !
  public:: Equil_FasTrace_EnTete
  public:: Equil_FasTrace_Write
  
contains

subroutine Equil_Errors_Warning(iErr)
  integer,intent(in):: iErr
  !
  if (iErr>=0) then
  
    call Info_Value("Error Code", iErr)
    call Warning_("Equil_Errors Called with iErr >= 0 !!")
  
  else
    
    select case(iErr)
    case(-1)  ;  call Warning_("Equil_Error: NewtMaxIts Exceeded in Newton")
    case(-2)  ;  call Warning_("Equil_Error: Singular Matrix in Newton")
    case(-3)  ;  call Warning_("Equil_Error: Roundoff Problem in LineSearch")
    case(-4)  ;  call Warning_("Equil_Error: No Convergence on Gammas")
    case(-5)  ;  call Warning_("Equil_Error: Stationary Point in Newton")
    case(-7)  ;  call Warning_("Equil_Error: MaxIter Exceeded in Outer Loop")
    case default
      call Info_Value("Error Code", iErr)
      call Warning_("Equil_Error: Unknown Code")
    end select
    
  end if
end subroutine Equil_Errors_Warning

subroutine Equil_Errors(iErr)
  integer,intent(in):: iErr
  !
  if (iErr>=0) then
  
    call Info_Value("Error Code", iErr)
    call Warning_("Equil_Errors Called with iErr >= 0 !!")
  
  else
    
    select case(iErr)
    case(-1)  ;  call Stop_("Equil_Error: NewtMaxIts Exceeded in Newton")
    case(-2)  ;  call Stop_("Equil_Error: Singular Matrix in Newton")
    case(-3)  ;  call Stop_("Equil_Error: Roundoff Problem in LineSearch")
    case(-4)  ;  call Stop_("Equil_Error: No Convergence on Gammas")
    case(-5)  ;  call Stop_("Equil_Error: Stationary Point in Newton")
    case(-7)  ;  call Stop_("Equil_Error: MaxIter Exceeded in Outer Loop")
    case default
      call Info_Value("Error Code", iErr)
      call Stop_("Equil_Error: Unknown Code")
    end select
    
  end if
end subroutine Equil_Errors

subroutine Equil_Zero(Cod)
!--
!-- for "very initial" speciation
!--
  use M_SolModel_Tools,only: Solmodel_Init_TotO_TotH
  use M_Equil_Read ! <-equil_Read_Debug,Equil_Read_Conditions,Equil_Read_PhaseAdd,Equil_Read_YesList
  use M_Equil_Vars, only: Equil_Vars_Alloc
  !
  use M_Global_Vars,only: vSpc,vEle,vFas,nAq,vMixModel
  use M_Global_Vars,only: SolModel
  use M_System_Vars,only: vCpn
  use M_Basis_Vars, only: iO_,iH_,iBal,isW,MWSv,initMolNu,bOH2
  use M_Basis_Vars, only: vLAx,vLCi,nAs,tAlfFs
  use M_Equil_Vars, only: DebFormula,DebNewt,DebJacob,vYesList
  use M_Equil_Vars, only: LogForAqu,DirectSub,cMethod,bFinDif
  use M_Equil_Vars, only: vLnAct
  use M_Equil_Vars, only: vTotS
  use M_Equil_Vars, only: nEquFas
  !
  character(len=3),intent(in):: Cod
  !!type(T_Component),dimension(:),intent(inout):: vCpn
  !
  integer,parameter:: AdjustMode= 1
  !
  logical:: Error
  character(len=80):: ErrorMsg
  integer:: I,iSpc
  !
  logical:: test_initmolnu= .false.
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Equil_Zero"
  !
  !--------------------------------------------------- default values --
  bFinDif=   .false.
  LogForAqu= .true.
  DirectSub= .false.  !! .true.  !! 
  bOH2=      .false.
  initMolNu= 1.0D-6
  !--------------------------------------------------/ default values --
  !
  !-------------------------------------------------- read CONDITIONS --
  call Equil_Read_Debug(DebFormula,DebNewt,DebJacob) !OUT
  call Equil_Read_Conditions(initMolNu,bOH2)  ! to be replaced by Equil_Read_Numeric !!
  call Equil_Read_Numeric( &
  & initMolNu,bOH2,LogForAqu,DirectSub,cMethod,bFinDif,Error,ErrorMsg)
  !
  !if(.not. LogForAqu) DirectSub= .true.
  !-------------------------------------------------/ read CONDITIONS --
  !
  call Solmodel_Init_TotO_TotH( & !
  & AdjustMode,        &          !IN
  & iBal,iO_,iH_,vEle, &          !IN
  & vCpn)                         !INOUT
  !-> calc. mole.nrs of O and H in fluid from mole.nrs of other elements
  !
  vFas(:)%MolFs= Zero
  nEquFas= 0
  !
  if(Cod(1:2)=="EQ") call Equil_Read_PhaseAdd(vFas)
  !
  !-> read SYSTEM.ROCK -> modify vFas(:)%Mole
  !~ ! new deleted 18/10/2007 18:25
  !~ do I=1,size(vFas)
    !~ vCpn(:)%Mole= vCpn(:)%Mole + tAlfFs(:,I) *vFas(I)%Mole
  !~ end do
  !
  !------------------------------------------------------ allocations --
  call Equil_Vars_Alloc(size(vEle),size(vSpc),nAq,nAs,size(vFas))
  !! print *,"vLnAct(isW)=", vLnAct(isW)  ;  call pause_
  !-----------------------------------------------------/ allocations --
  !
  !! if(trim(Cod)/="SPC") call Equil_Read_YesList(vFas,vYesList)
  if(Cod(1:2)=="EQ") call Equil_Read_YesList(size(vFas),vFas,vYesList)
  !
  !-------------------------------------------------- initial gamma's --
  vSpc(:)%Dat%LGam= Zero
  !
  !------------------------------------ initial mole numbers in fluid --
  != -> initialize vSpc(1:nAq)%Mole, used as initial guess for solver
  !Equil_ZeroMolF !-> initial vSpc(1:nAq)%Dat%Mole
  vSpc(1:nAq)%Dat%Mole= InitMolNu !initial value for all aqu'species
  vSpc(isW)%Dat%Mole=   One/MWSv  !except Solvent
  !
  if(test_initmolnu) then
    do I=1,size(vCpn)
      if(I/=iH_) then
        iSpc= vCpn(I)%iSpc
        if(vLCi(iSpc)) vSpc(iSpc)%Dat%Mole= vCpn(I)%Mole
      end if
    end do
  end if
  !
  != for mobile aqu'species, compute mole numbers estimates
  != from their activities (which are fixed) and the current gamma's
  do I=1,nAq
    if(vLAx(I)) &
    & vSpc(I)%Dat%Mole= exp(vSpc(I)%Dat%LAct -vSpc(I)%Dat%LGam)
    !this gives molality, should multiply by mass of solvent ...
  end do
  !-----------------------------------/ initial mole numbers in fluid --
  !!
  !call Specia_InitBalance(vTotS) !!NEW201010!!
  if(iBal/=0) vTotS(iBal)= Zero
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Equil_Zero"
  !
end subroutine Equil_Zero

subroutine Equil_Restart
!--
!-- restart a speciation, using vCpn%Mole, vSpc%Dat%LnAct, vSpc%LnGam
!--
  use M_IoTools,    only: OutStrVec
  !
  use M_Global_Vars,only: vSpc,vFas,nAq
  
  use M_System_Vars,only: vCpn
  
  use M_Basis_Vars, only: iBal,nCi,nAs,tNuAs,tNuFas,tAlfFs,vOrdPr,vOrdAs
  
  use M_Equil_Vars, only: vDeltaG_As,vTotS,vLnAct,vLnGam,vMolF
  use M_Equil_Vars, only: vDeltaG_Fs,vFasMole,vAffScale
  !
  integer:: I
  !
  !! DEBUGG !!
  ! write(6,'(A)') "==Equil_Restart=="
  ! do I=1,size(vFas)
  !   write(6,'(A)')  trim(vFas(I)%NamFs)
  ! end do
  ! pause

  !------------------------------------------------- update deltaG's  --
  ! update deltaG's of second'aqu'species (may have changed if TP changed ...)
  do I=1,nAs
    vDeltaG_As(I)= &
    & vSpc(vOrdAs(I))%G0rt - dot_product(tNuAs(I,:),vSpc(vOrdPr(:))%G0rt)
  end do
  !
  ! update deltaG's of phases (may have changed if TP changed ...)
  ! DeltaG is free energy change for G(phase) - G(prim'species)
  ! i.e. DeltaG of reaction of formatION of phase from prim'species
  do I=1,size(vFas)
    vDeltaG_Fs(I)= &
    & vFas(I)%Grt - dot_product(tNuFas(I,:),vSpc(vOrdPr(:))%G0rt)
  end do
  !------------------------------------------------/ update deltaG's  --
  !
  do I=1,size(vFas)
    vAffScale(I)= SUM(ABS(tNuFas(I,:)))
  end do
  !
  ! retrieve data from current vSpc 
  vLnAct(:)=    vSpc(:)%Dat%LAct
  vLnGam(:)=    vSpc(:)%Dat%LGam
  vMolF(1:nAq)= vSpc(1:nAq)%Dat%Mole
  !
  !~ print *,"Equil_Restart"
  !~ do i=1,size(vFas)
    !~ if(vFas(i)%Mole>0.D0) &
    !~ & print *,"=========",vFas(i)%Mole,"=",trim(vFas(i)%NamFs)
  !~ end do
  !~ pause
  
  ! retrieve data from current vFas
  vFasMole(:)=  vFas(:)%MolFs
  !
  ! vTotS= mole nr components in SYSTEM (fluid + min), RHS in residual
  ! == in SPC, vTotS = vTot_inFluid
  ! == in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
  !
  ! retrieve data from current vCpn
  vTotS(1:nCi)= vCpn(1:nCi)%Mole ! mole nrs in fluid
  !---------------------------------------------------------- REWORK !!!
  != add mole nrs from non-fluid (case of PATHEQ* calculations)
  do I=1,size(vFas)
    vTotS(1:nCi)= vTotS(1:nCi) + tAlfFs(1:nCi,I) *vFasMole(I)
  end do
  !! if(iDebug>2) call OutStrVec(51,vTotS(1:nCi))
  
  if(iBal/=0) vTotS(iBal)= Zero
  
  !print *,"iBal=",iBal   ;  pause
  !call Specia_InitBalance(vTotS) !!NEW201010!!
  
  return
end subroutine Equil_Restart

subroutine Specia_InitBalance(vTot_)
  use M_Basis,     only: Basis_InitBalance
  use M_Basis_Vars,only: iBal,vOrdPr,vOrdAq,vOrdAs,vOrdMs
  use M_Basis_Vars,only: tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs
  !
  real(dp),dimension(:),intent(inout):: vTot_
  !
  !! if H is "INERT", mole balance on H replaced by mole balance on charges  !
  !! if(trim(vCpn(iH_)%Statut)=="INERT") iBal=iH_ !NEW_17/10/2007 17:40
  !
  if(iBal/=0) call Basis_InitBalance( & !
  & iBal,                 & !IN
  & vOrdPr,vOrdAs,vOrdMs, & !IN
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !INOUT
  & tAlfFs)                 !INOUT
  !
  if(iBal/=0) vTot_(iBal)= Zero
  !
end subroutine Specia_InitBalance

subroutine Equil_Save(lEquil)
!--
!-- Update vSpc & vCpn for use in next speciations
!--
  use M_Global_Vars,only: vEle,vSpc,vFas,nAq
  use M_System_Vars,only: vCpn
  use M_Basis_Vars, only: nCi,nAx,nAs,tAlfPr,tAlfAs,tAlfFs,vOrdPr,vOrdAs
  use M_Equil_Vars, only: vMolF,vLnAct,vLnGam,vFasMole,nEquFas
  !
  use M_Numeric_Const,only: Ln10
  use M_Basis_Vars,   only: isH_,isOH
  !
  logical,intent(in),optional:: lEquil
  !
  logical :: B
  integer :: iCp
  real(dp):: R
  
  B=.false.
  if(present(lEquil)) B=lEquil
  
  !-------------- store results in vSpc to use in subsequent calcul's --
  vSpc(1:nAq)%Dat%Mole= vMolF(1:nAq) !/vMolF(1)/MWsv
  vSpc(:)%Dat%LAct=     vLnAct(:) !_____id
  vSpc(:)%Dat%LGam=     vLnGam(:) !_____id
  
  !-------------- store results in vFas to use in subsequent calcul's --
  !~ vFas(:)%Mole= vFasMole(:)
  if(nEquFas>0) call Save_EquPhase
  
  !do iCp=1,size(vCpn) !not useful here ??
  !  if(vCpn(iCp)%Statut/="INERT") vCpn(iCp)%Mole= & 
  !  & SUM(tAlfPr(iCp,1    :nCi    ) *vMolF(vOrdPr(1:nCi)))  & !"internal" prim'species
  !  + SUM(tAlfPr(iCp,nCi+1:nCi+nAx) *vMolF(vOrdPr(nCi+1:nCi+nAx))) &
  !  + SUM(tAlfAs(iCp,1    :nAs    ) *vMolF(vOrdAs(1:nAs)))
  !end do
  
  !------------- compute mole nrs components in fluid -> vCpn(:)%Mole --
  do iCp=1,size(vCpn) !
    !
    R= SUM(tAlfPr(iCp,1    :nCi    ) *vMolF(vOrdPr(1:nCi)))  &
    +  SUM(tAlfPr(iCp,nCi+1:nCi+nAx) *vMolF(vOrdPr(nCi+1:nCi+nAx))) &
    +  SUM(tAlfAs(iCp,1    :nAs    ) *vMolF(vOrdAs(1:nAs)))
    !
    vCpn(iCp)%Mole= R
    !
  end do
  
  if(iDebug==4) then
    print '(A)',"Equil_Save, check components' mole nrs"
    do iCp=1,size(vCpn) !
      print '(2A,2G15.6)',&
      & "Equil_Save, ", vEle(vCpn(iCp)%iEle)%NamEl,vCpn(iCp)%Mole
    end do
    call Pause_
  end if
  
  return
end subroutine Equil_Save

subroutine Save_EquPhase
  use M_T_MixModel, only: MaxPole
  use M_T_Phase,    only: T_Phase
  use M_Equil_Vars, only: T_EquPhase
  use M_Global_Vars,only: vMixModel,vMixFas,vFas
  use M_Basis_Vars, only: tAlfFs,tNuFas
  use M_Equil_Vars, only: nEquFas,vEquFas,cEquMode
  !
  type(T_Phase),allocatable:: vFas0(:)
  real(dp),     allocatable:: tAlf(:,:),tNu(:,:)
  integer:: vIPole(MaxPole) ! indexes of end-members in vFas
  integer:: nPur,nMix
  integer:: I,K,nP,C,nCp
  
  nPur= count(vFas(:)%iSpc>0)
  nCp= size(tAlfFs,1)
  
  !--- build vMixFas --
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
        nP= vMixModel(vEquFas(I)%iMix)%nPole
        vMixFas(K)%Name= trim(vEquFas(I)%NamEq)
        vMixFas(K)%iModel= 1 !vEquFas(I)%iMix
        vMixFas(K)%vXPole(1:nP)= vEquFas(I)%vXPole(1:nP)
        
      end if
    end do
  end if
  !--- build vMixFas --
  
  nMix= count(vEquFas(1:nEquFas)%iMix>0)
  
  if(nMix>0) then
  
    ! temporary copy
    allocate(vFas0(nPur))
    vFas0(1:nPur)= vFas(1:nPur)
    
    deallocate(vFas)
    allocate(vFas(nPur+nMix))
    vFas(1:nPur)= vFas0(1:nPur)
    
    deallocate(vFas0)
  
    allocate(tAlf(nCp,nPur))        ; tAlf(:,1:nPur)= tAlfFs(:,1:nPur)
    deallocate(tAlfFs)
    allocate(tAlfFs(nCp,nPur+nMix)) ; tAlfFs(:,1:nPur)= tAlf(:,1:nPur)
    deallocate(tAlf)

    allocate(tNu(nPur,nCp))         ; tNu(1:nPur,:)= tNuFas(1:nPur,:)
    deallocate(tNuFas)
    allocate(tNuFas(nPur+nMix,nCp)) ; tNuFas(1:nPur,:)= tNu(1:nPur,:)
    deallocate(tNu)

  end if
  
  vFas(:)%MolFs= Zero
  
  K=nPur
  do I=1,nEquFas
  
    if(vEquFas(I)%iPur>0) vFas(vEquFas(I)%iPur)%MolFs= vEquFas(I)%MolFs
    
    if(vEquFas(I)%iMix>0) then
    
      K= K+1
      
      vFas(K)%MolFs= vEquFas(I)%MolFs
      vFas(K)%iSpc=  0
      vFas(K)%iSol=  0
      vFas(K)%iMix=  K -nPur
      vFas(K)%NamFs= trim(vEquFas(I)%NamEq)
      !! vFas(K)%Typ=   "MIXT"
      
      nP= vEquFas(I)%NPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      do C=1,nCp
        tAlfFs(C,K)= SUM( vEquFas(I)%vXPole(1:nP) &
        &          * tAlfFs(C,vIPole(1:nP)) )
        tNuFas(K,C)= SUM( vEquFas(I)%vXPole(1:nP) &
        &          * tNuFas(vIPole(1:nP),C) )
      end do

    end if
    
  end do
  
  return
end subroutine Save_EquPhase

subroutine Equil_CalcFasAff
!--
!-- calculate affinity of pure phases (i.e. non-aqu'species)
!--
  use M_Global_Vars,only: vSpc,vFas
  use M_Basis_Vars, only: tNuFas,vOrdPr
  use M_Equil_Vars, only: vFasAff
  !
  integer  :: iFs,nFs,nCp
  
  nFs= size(vFas)
  nCp= size(vOrdPr)
  
  if (nFs>0) then
    !note:
    !vDeltaG_Ms(iMs)=  vSpc(vOrdMs(iMs))%G0rt 
    !               - dot_product(tNuMs(iMs,1:nCp),Spc(vOrdPr(1:nCp))%G0rt) 
    !
    !-> deltaG of PRECIPITATION reaction, Prim'Species -> Mineral
    !
    !affinity/RT of formatION reaction of mineral,
    != A/RT= -DotProd(Nu,Mu)= -deltaG of Prim'Sp->Mineral, A/RT=ln(K/Q)
    !
    !LnQsK(iMs)= dot_product( tNuMs(iMs,1:nCp), vSpc(vOrdPr(1:nCp))%LnAct ) 
    !          - vDeltaG_Ms(iMs)
    !
    do iFs=1,nFs
      vFasAff(iFs)= &
      & vFas(iFs)%Grt &
      - dot_product( tNuFas(iFs,1:nCp), &
      &              vSpc(vOrdPr(1:nCp))%Dat%LAct+vSpc(vOrdPr(1:nCp))%G0rt )
      !vFasAff(iFs)= vFasAff(iFs) /SUM(ABS(tNuFas(iFs,1:nCp)))
    end do !_______________________________end
  end if
  
  !expression used in Equil_1:
  !
  !  vFasAff(iFs)= vDeltaG_Fs(iFs) &
  !  &           - dot_product(tNuFas(iFs,1:nCp),vLnAct(vOrdPr(1:nCp)))
  !
  !with
  !
  !  vDeltaG_Fs(:)= vFas(:)%Grt 
  !               - matmul(tNuFas(:,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
  !
  !-> can also calculate vFasAff as
  !
  !  vFasAff(iFs)= &
  !  & vSpc(vFas(iFs)%iSpc)%G0rt &
  !  - dot_product(tNuFas(iFs,1:nCp), &
  !  &             vSpc(vOrdPr(1:nCp))%Dat%LAct+vSpc(vOrdPr(1:nCp))%G0rt)
  
  return
end subroutine Equil_CalcFasAff

subroutine Equil_Trace_Close
!--
!-- close special trace files for activities, Newton, etc
!--
  use M_Equil_Vars, only: fTrAct, fActiz, fTrcEq
  use M_Numeric_Tools,only: fNewtF,fNewtR
  !
  if(fTrAct>0)  then; close(fTrAct); fTrAct= 0; end if
  if(fActiZ>0)  then; close(fActiZ); fActiZ= 0; end if
  if(fTrcEq>0)  then; close(fTrcEq); fTrcEq= 0; end if
  !
  if(fNewtF>0)  then; close(fNewtF); fNewtF= 0; end if
  if(fNewtR>0)  then; close(fNewtR); fNewtR= 0; end if
  !
  return
end subroutine Equil_Trace_Close

subroutine Equil_Trace_Init
!--
!-- open special trace files for activities, Newton, etc
!--
  use M_IOTools,    only: GetUnit
  use M_Files,      only: DirOut,Files_Index_Write
  use M_Numeric_Tools,only: fNewtF,fNewtR,fNewt_I
  !
  use M_Global_Vars,only: vSpc
  use M_Equil_Write,only: Equil_Write_EnTete
  use M_Basis_Vars, only: vLCi,vLAs,vLAx,vPrmFw
  use M_Equil_Vars, only: DebActiv,fTrAct,fActiZ,DebNewt
  !
  integer:: I
  
  if(DebActiv) then
    !
    call GetUnit(fTrAct)
    open(fTrAct,file=trim(DirOut)//"_activ_deb.log")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_activ_deb.log",&
    & "EQUIL/LOG: activ' aqu'species")
    write(fTrAct,'(A,A1)',advance="no") "indx",T_
    do I=1,size(vSpc)
      if(vSpc(I)%Typ=="AQU") write(fTrAct,'(A,A1)',advance="no") trim(vSpc(I)%NamSp),T_
    end do
    write(fTrAct,*)
    !
    call GetUnit(fActiz); open(fActiz,file=trim(DirOut)//"_activ.log")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_activ.log",&
    & "EQUIL/LOG: log10(Activ)")
    call Equil_Write_EnTete(fActiz,vSpc,vPrmFw,vLCi.or.vLAs.or.vLAx,Str1="indx")
    !
  else
    !
    fTrAct= 0
    fActiZ= 0
    !
  end if
  
  if(DebNewt) then
    !
    call GetUnit(fNewtF)
    open(fNewtF,file=trim(DirOut)//"_newtequspc1.log")
    call Equil_Write_EnTete(fNewtF,vSpc,vPrmFw,vLCi.or.vLAs,Str1="count1",Str2="count2")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_newtequspc1.log",&
    & "EQUIL/LOG: Newton solution")
    !
    call GetUnit(fNewtR)
    open(fNewtR,file=trim(DirOut)//"_newtequspc2.log")
    call Equil_Write_EnTete(fNewtR,vSpc,vPrmFw,vLCi.or.vLAs,Str1="count1",Str2="count2")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_newtequspc2.log",&
    & "EQUIL/LOG: Newton residue")
    !
    fNewt_I=0
    !
  else
    !
    fNewtF=0
    fNewtR=0
    !
  end if
  
  return
end subroutine Equil_Trace_Init

subroutine Equil_FasTrace_EnTete(vFas,vYesList,F)

  use M_IOTools,only: GetUnit
  use M_Files,  only: DirOut
  use M_T_Phase,only: T_Phase
  !
  type(T_Phase),intent(in) :: vFas(:)
  logical,      intent(in) :: vYesList(:)
  integer,      intent(out):: F
  !
  integer:: iFs,nFs
  
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_equil.log")
  !F=51
  
  write(F,'(3(A,A1))',advance="NO") "iCount",T_,"nIts",T_,"dXi",T_
  
  nFs=size(vFas)
  !
  do iFs=1,nFs
    if(vYesList(iFs)) &
    & write(F,'(A,A1)',advance="NO") "Mol_"//trim(vFas(iFs)%NamFs),T_
  end do
  do iFs=1,nFs
    if(vYesList(iFs)) &
    & write(F,'(A,A1)',advance="NO") "lQsK_"//trim(vFas(iFs)%NamFs),T_
  end do
  write(F,*)
  
  return
end subroutine Equil_FasTrace_EnTete

subroutine Equil_FasTrace_Write(F,ICount,nIts,dXi,vYesList,vFasMole,vFasAff)
  use M_Numeric_Const,  only: Ln10
  use M_IOTools,only: OutStrVec
  !
  integer, intent(in):: F
  integer, intent(in):: iCount
  integer, intent(in):: nIts
  real(dp),intent(in):: dXi
  logical, dimension(:),intent(in):: vYesList
  real(dp),dimension(:),intent(in):: vFasMole,vFasAff
  
  write(F,'(2(I7,A1),G15.6,A1)',advance="NO") iCount,T_,nIts,T_,dXi,T_
  call OutStrVec(F,vFasMole,     vL=vYesList,CR=.false.)
  call OutStrVec(F,-vFasAff/Ln10,vL=vYesList)
  
  return
end subroutine Equil_FasTrace_Write

end module M_Equil_Tools

