MODULE M_Equil_Tools
!--
!-- common tools for equilibrium calculations
!--
  USE M_Trace,ONLY: iDebug,fTrc,fHtm,T_,Stop_,Pause_,Warning_
  USE M_Info_Value
  USE M_Kinds
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  !!PUBLIC:: Equil_Solve
  !
  PUBLIC:: Equil_Zero      !for a new speciation or equilibrium
  PUBLIC:: Equil_Restart
  PUBLIC:: Equil_Save      !store results of Equil_Calc in vCpn, vSpc 
  PUBLIC:: Equil_Errors
  PUBLIC:: Equil_Errors_Warning
  !
  PUBLIC:: Equil_CalcFasAff
  !
  PUBLIC:: Equil_Trace_Init
  PUBLIC:: Equil_Trace_Close
  !
  PUBLIC:: Equil_FasTrace_EnTete
  PUBLIC:: Equil_FasTrace_Write
  
CONTAINS

SUBROUTINE Equil_Errors_Warning(iErr)
  INTEGER,INTENT(IN):: iErr
  !
  IF (iErr>=0) THEN
  
    CALL Info_Value("Error Code", iErr)
    CALL Warning_("Equil_Errors Called with iErr >= 0 !!")
  
  ELSE
    
    SELECT CASE(iErr)
    CASE(-1)  ;  CALL Warning_("Equil_Error: NewtMaxIts Exceeded in Newton")
    CASE(-2)  ;  CALL Warning_("Equil_Error: Singular Matrix in Newton")
    CASE(-3)  ;  CALL Warning_("Equil_Error: Roundoff Problem in LineSearch")
    CASE(-4)  ;  CALL Warning_("Equil_Error: No Convergence on Gammas")
    CASE(-5)  ;  CALL Warning_("Equil_Error: Stationary Point in Newton")
    CASE(-7)  ;  CALL Warning_("Equil_Error: MaxIter Exceeded in Outer Loop")
    CASE DEFAULT
      CALL Info_Value("Error Code", iErr)
      CALL Warning_("Equil_Error: Unknown Code")
    ENDSELECT
    
  END IF
ENDSUBROUTINE Equil_Errors_Warning

SUBROUTINE Equil_Errors(iErr)
  INTEGER,INTENT(IN):: iErr
  !
  IF (iErr>=0) THEN
  
    CALL Info_Value("Error Code", iErr)
    CALL Warning_("Equil_Errors Called with iErr >= 0 !!")
  
  ELSE
    
    SELECT CASE(iErr)
    CASE(-1)  ;  CALL Stop_("Equil_Error: NewtMaxIts Exceeded in Newton")
    CASE(-2)  ;  CALL Stop_("Equil_Error: Singular Matrix in Newton")
    CASE(-3)  ;  CALL Stop_("Equil_Error: Roundoff Problem in LineSearch")
    CASE(-4)  ;  CALL Stop_("Equil_Error: No Convergence on Gammas")
    CASE(-5)  ;  CALL Stop_("Equil_Error: Stationary Point in Newton")
    CASE(-7)  ;  CALL Stop_("Equil_Error: MaxIter Exceeded in Outer Loop")
    CASE DEFAULT
      CALL Info_Value("Error Code", iErr)
      CALL Stop_("Equil_Error: Unknown Code")
    ENDSELECT
    
  END IF
ENDSUBROUTINE Equil_Errors

SUBROUTINE Equil_Zero(Cod)
!--
!-- for "very initial" speciation
!--
  USE M_SolModel_Tools,ONLY: Solmodel_Init_TotO_TotH
  USE M_Equil_Read ! <-equil_Read_Debug,Equil_Read_Conditions,Equil_Read_PhaseAdd,Equil_Read_YesList
  USE M_Equil_Vars, ONLY: Equil_Vars_Alloc
  !
  USE M_Global_Vars,ONLY: vSpc,vEle,vFas,nAq,vMixModel
  USE M_Global_Vars,ONLY: SolModel
  USE M_System_Vars,ONLY: vCpn
  USE M_Basis_Vars, ONLY: iO_,iH_,iBal,isW,MWSv,initMolNu,bOH2
  USE M_Basis_Vars, ONLY: vLAx,vLCi,nAs,tAlfFs
  USE M_Equil_Vars, ONLY: DebFormula,DebNewt,DebJacob,vYesList
  USE M_Equil_Vars, ONLY: LogForAqu,DirectSub,cMethod,bFinDif
  USE M_Equil_Vars, ONLY: vTotS
  USE M_Equil_Vars, ONLY: nEquFas
  !
  CHARACTER(LEN=3),INTENT(IN):: Cod
  !!TYPE(T_Component),DIMENSION(:),INTENT(INOUT):: vCpn
  !
  INTEGER,PARAMETER:: AdjustMode= 1
  !
  LOGICAL:: Error
  CHARACTER(LEN=80):: ErrorMsg
  INTEGER:: I,iSpc
  !
  LOGICAL:: test_initmolnu= .FALSE.
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Zero"
  !
  !--------------------------------------------------- default values --
  bFinDif=   .FALSE.
  LogForAqu= .TRUE.
  DirectSub= .FALSE.  !! .TRUE.  !! 
  bOH2=      .FALSE.
  initMolNu= 1.0D-6
  !--------------------------------------------------/ default values --
  !
  !-------------------------------------------------- read CONDITIONS --
  CALL Equil_Read_Debug(DebFormula,DebNewt,DebJacob) !OUT
  CALL Equil_Read_Conditions(initMolNu,bOH2)  ! to be replaced by Equil_Read_Numeric !!
  CALL Equil_Read_Numeric( &
  & initMolNu,bOH2,LogForAqu,DirectSub,cMethod,bFinDif,Error,ErrorMsg)
  !
  !IF(.NOT. LogForAqu) DirectSub= .TRUE.
  !-------------------------------------------------/ read CONDITIONS --
  !
  CALL Solmodel_Init_TotO_TotH( & !
  & AdjustMode,        &          !IN
  & iBal,iO_,iH_,vEle, &          !IN
  & vCpn)                         !INOUT
  !-> calc. mole.nrs of O and H in fluid from mole.nrs of other elements
  !
  vFas(:)%Mole=Zero
  nEquFas= 0
  !
  IF(Cod(1:2)=="EQ") CALL Equil_Read_PhaseAdd(vFas)
  !
  !-> READ SYSTEM.ROCK -> modify vFas(:)%Mole
  !~ ! new deleted 18/10/2007 18:25
  !~ DO I=1,SIZE(vFas)
    !~ vCpn(:)%Mole= vCpn(:)%Mole + tAlfFs(:,I) *vFas(I)%Mole
  !~ ENDDO
  !
  !------------------------------------------------------ allocations --
  CALL Equil_Vars_Alloc(SIZE(vEle),SIZE(vSpc),nAq,nAs,SIZE(vFas))
  !-----------------------------------------------------/ allocations --
  !
  !! IF(TRIM(Cod)/="SPC") CALL Equil_Read_YesList(vFas,vYesList)
  IF(Cod(1:2)=="EQ") CALL Equil_Read_YesList(SIZE(vFas),vFas,vYesList)
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
  IF(test_initmolnu) THEN
    DO I=1,SIZE(vCpn)
      IF(I/=iH_) THEN
        iSpc= vCpn(I)%iSpc
        IF(vLCi(iSpc)) vSpc(iSpc)%Dat%Mole= vCpn(I)%Mole
      ENDIF
    ENDDO
  ENDIF
  !
  != for mobile aqu'species, compute mole numbers estimates
  != from their activities (which are fixed) and the current gamma's
  DO I=1,nAq
    IF(vLAx(I)) &
    & vSpc(I)%Dat%Mole= EXP(vSpc(I)%Dat%LAct -vSpc(I)%Dat%LGam)
    !this gives molality, should multiply by mass of solvent ...
  ENDDO
  !-----------------------------------/ initial mole numbers in fluid --
  !!
  !CALL Specia_InitBalance(vTotS) !!NEW201010!!
  IF(iBal/=0) vTotS(iBal)= Zero
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Zero"
  !
ENDSUBROUTINE Equil_Zero

SUBROUTINE Equil_Restart
!--
!-- restart a speciation, using vCpn%Mole, vSpc%Dat%LnAct, vSpc%LnGam
!--
  USE M_IoTools,    ONLY: OutStrVec
  !
  USE M_Global_Vars,ONLY: vSpc,vFas,nAq
  
  USE M_System_Vars,ONLY: vCpn
  
  USE M_Basis_Vars, ONLY: iBal,nCi,nAs,tNuAs,tNuFas,tAlfFs,vOrdPr,vOrdAs
  
  USE M_Equil_Vars, ONLY: vDeltaG_As,vTotS,vLnAct,vLnGam,vMolF
  USE M_Equil_Vars, ONLY: vDeltaG_Fs,vFasMole,vAffScale
  !
  INTEGER:: I
  !
  !! DEBUGG !!
  ! WRITE(6,'(A)') "==Equil_Restart=="
  ! DO I=1,SIZE(vFas)
  !   WRITE(6,'(A)')  TRIM(vFas(I)%NamFs)
  ! END DO
  ! PAUSE

  !------------------------------------------------- update deltaG's  --
  ! update deltaG's of second'aqu'species (may have changed if TP changed ...)
  DO I=1,nAs
    vDeltaG_As(I)= &
    & vSpc(vOrdAs(I))%G0rt - DOT_PRODUCT(tNuAs(I,:),vSpc(vOrdPr(:))%G0rt)
  ENDDO
  !
  ! update deltaG's of phases (may have changed if TP changed ...)
  ! DeltaG is free energy change for G(phase) - G(prim'species)
  ! i.e. DeltaG of reaction of FORMATION of phase from prim'species
  DO I=1,SIZE(vFas)
    vDeltaG_Fs(I)= &
    & vFas(I)%Grt - DOT_PRODUCT(tNuFas(I,:),vSpc(vOrdPr(:))%G0rt)
  ENDDO
  !------------------------------------------------/ update deltaG's  --
  !
  DO I=1,SIZE(vFas)
    vAffScale(I)= SUM(ABS(tNuFas(I,:)))
  ENDDO
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
  !~ enddo
  !~ pause
  
  ! retrieve data from current vFas
  vFasMole(:)=  vFas(:)%Mole
  !
  ! vTotS= mole nr components in SYSTEM (fluid + min), RHS in residual
  ! == in SPC, vTotS = vTot_inFluid
  ! == in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
  !
  ! retrieve data from current vCpn
  vTotS(1:nCi)= vCpn(1:nCi)%Mole ! mole nrs in fluid
  !---------------------------------------------------------- REWORK !!!
  != add mole nrs from non-fluid (case of PATHEQ* calculations)
  DO I=1,SIZE(vFas)
    vTotS(1:nCi)= vTotS(1:nCi) + tAlfFs(1:nCi,I) *vFas(I)%Mole
  ENDDO
  !! IF(iDebug>2) CALL OutStrVec(51,vTotS(1:nCi))
  
  IF(iBal/=0) vTotS(iBal)= Zero
  
  !print *,"iBal=",iBal   ;  pause
  !CALL Specia_InitBalance(vTotS) !!NEW201010!!
  
  RETURN
ENDSUBROUTINE Equil_Restart

SUBROUTINE Specia_InitBalance(vTot_)
  USE M_Basis,     ONLY: Basis_InitBalance
  USE M_Basis_Vars,ONLY: iBal,vOrdPr,vOrdAq,vOrdAs,vOrdMs
  USE M_Basis_Vars,ONLY: tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs
  !
  REAL(dp),DIMENSION(:),INTENT(INOUT):: vTot_
  !
  !! IF H is "INERT", mole balance on H replaced by mole balance on charges  !
  !! IF(TRIM(vCpn(iH_)%Statut)=="INERT") iBal=iH_ !NEW_17/10/2007 17:40
  !
  IF(iBal/=0) CALL Basis_InitBalance( & !
  & iBal,                 & !IN
  & vOrdPr,vOrdAs,vOrdMs, & !IN
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !INOUT
  & tAlfFs)                 !INOUT
  !
  IF(iBal/=0) vTot_(iBal)= Zero
  !
ENDSUBROUTINE Specia_InitBalance

SUBROUTINE Equil_Save(lEquil)
!--
!-- Update vSpc & vCpn for use in next speciations
!--
  USE M_Global_Vars,ONLY: vEle,vSpc,vFas,nAq
  USE M_System_Vars,ONLY: vCpn
  USE M_Basis_Vars, ONLY: nCi,nAx,nAs,tAlfPr,tAlfAs,tAlfFs,vOrdPr,vOrdAs
  USE M_Equil_Vars, ONLY: vMolF,vLnAct,vLnGam,vFasMole,nEquFas
  !
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Basis_Vars,   ONLY: isH_,isOH
  !
  LOGICAL,INTENT(IN),OPTIONAL:: lEquil
  !
  LOGICAL :: B
  INTEGER :: iCp
  REAL(dp):: R
  
  B=.FALSE.
  IF(PRESENT(lEquil)) B=lEquil
  
  !-------------- store results in vSpc to use in subsequent calcul's --
  vSpc(1:nAq)%Dat%Mole= vMolF(1:nAq) !/vMolF(1)/MWsv
  vSpc(:)%Dat%LAct=     vLnAct(:) !_____id
  vSpc(:)%Dat%LGam=     vLnGam(:) !_____id
  
  !-------------- store results in vFas to use in subsequent calcul's --
  !~ vFas(:)%Mole= vFasMole(:)
  IF(nEquFas>0) CALL Save_EquPhase
  
  !DO iCp=1,SIZE(vCpn) !not useful here ??
  !  IF(vCpn(iCp)%Statut/="INERT") vCpn(iCp)%Mole= & 
  !  & SUM(tAlfPr(iCp,1    :nCi    ) *vMolF(vOrdPr(1:nCi)))  & !"internal" prim'species
  !  + SUM(tAlfPr(iCp,nCi+1:nCi+nAx) *vMolF(vOrdPr(nCi+1:nCi+nAx))) &
  !  + SUM(tAlfAs(iCp,1    :nAs    ) *vMolF(vOrdAs(1:nAs)))
  !ENDDO
  
  !------------- compute mole nrs components in fluid -> vCpn(:)%Mole --
  DO iCp=1,SIZE(vCpn) !
    !
    R= SUM(tAlfPr(iCp,1    :nCi    ) *vMolF(vOrdPr(1:nCi)))  &
    +  SUM(tAlfPr(iCp,nCi+1:nCi+nAx) *vMolF(vOrdPr(nCi+1:nCi+nAx))) &
    +  SUM(tAlfAs(iCp,1    :nAs    ) *vMolF(vOrdAs(1:nAs)))
    !
    vCpn(iCp)%Mole= R
    !
  ENDDO
  
  IF(iDebug==4) THEN
    PRINT '(A)',"Equil_Save, check components' mole nrs"
    DO iCp=1,SIZE(vCpn) !
      PRINT '(2A,2G15.6)',&
      & "Equil_Save, ", vEle(vCpn(iCp)%iEle)%NamEl,vCpn(iCp)%Mole
    ENDDO
    CALL Pause_
  END IF
  
  RETURN
ENDSUBROUTINE Equil_Save

SUBROUTINE Save_EquPhase
  USE M_T_MixModel, ONLY: MaxPole
  USE M_T_Phase,    ONLY: T_Phase
  USE M_Equil_Vars, ONLY: T_EquPhase
  USE M_Global_Vars,ONLY: vMixModel,vMixFas,vFas
  USE M_Basis_Vars, ONLY: tAlfFs,tNuFas
  USE M_Equil_Vars, ONLY: nEquFas,vEquFas,cEquMode
  !
  TYPE(T_Phase),ALLOCATABLE:: vFas0(:)
  REAL(dp),     ALLOCATABLE:: tAlf(:,:),tNu(:,:)
  INTEGER:: vIPole(MaxPole) ! indexes of end-members in vFas
  INTEGER:: nPur,nMix
  INTEGER:: I,K,nP,C,nCp
  
  nPur= COUNT(vFas(:)%iSpc>0)
  nCp= SIZE(tAlfFs,1)
  
  !--- build vMixFas --
  K=0
  DO I=1,nEquFas
    IF(vEquFas(I)%iMix>0) K= K+1
  ENDDO
  DEALLOCATE(vMixFas)
  ALLOCATE(vMixFas(K))
  !
  K=0
  IF(SIZE(vMixFas)>0) THEN
    DO I=1,nEquFas
      IF(vEquFas(I)%iMix>0) THEN
      
        K= K+1
        nP= vMixModel(vEquFas(I)%iMix)%nPole
        vMixFas(K)%Name= TRIM(vEquFas(I)%NamEq)
        vMixFas(K)%iModel= 1 !vEquFas(I)%iMix
        vMixFas(K)%vXPole(1:nP)= vEquFas(I)%vXPole(1:nP)
        
      ENDIF
    ENDDO
  ENDIF
  !--- build vMixFas --
  
  nMix= COUNT(vEquFas(1:nEquFas)%iMix>0)
  
  IF(nMix>0) THEN
  
    ! temporary copy
    ALLOCATE(vFas0(nPur))
    vFas0(1:nPur)= vFas(1:nPur)
    
    DEALLOCATE(vFas)
    ALLOCATE(vFas(nPur+nMix))
    vFas(1:nPur)= vFas0(1:nPur)
    
    DEALLOCATE(vFas0)
  
    ALLOCATE(tAlf(nCp,nPur))        ; tAlf(:,1:nPur)= tAlfFs(:,1:nPur)
    DEALLOCATE(tAlfFs)
    ALLOCATE(tAlfFs(nCp,nPur+nMix)) ; tAlfFs(:,1:nPur)= tAlf(:,1:nPur)
    DEALLOCATE(tAlf)

    ALLOCATE(tNu(nPur,nCp))         ; tNu(1:nPur,:)= tNuFas(1:nPur,:)
    DEALLOCATE(tNuFas)
    ALLOCATE(tNuFas(nPur+nMix,nCp)) ; tNuFas(1:nPur,:)= tNu(1:nPur,:)
    DEALLOCATE(tNu)

  ENDIF
  
  vFas(:)%Mole= Zero
  
  K=nPur
  DO I=1,nEquFas
  
    IF(vEquFas(I)%iPur>0) vFas(vEquFas(I)%iPur)%Mole= vEquFas(I)%Mole
    
    IF(vEquFas(I)%iMix>0) THEN
    
      K= K+1
      
      vFas(K)%Mole= vEquFas(I)%Mole
      vFas(K)%iSpc=  0
      vFas(K)%iSol=  0
      vFas(K)%iMix=  K -nPur
      vFas(K)%NamFs= TRIM(vEquFas(I)%NamEq)
      !! vFas(K)%Typ=   "MIXT"
      
      nP= vEquFas(I)%NPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      DO C=1,nCp
        tAlfFs(C,K)= SUM( vEquFas(I)%vXPole(1:nP) &
        &          * tAlfFs(C,vIPole(1:nP)) )
        tNuFas(K,C)= SUM( vEquFas(I)%vXPole(1:nP) &
        &          * tNuFas(vIPole(1:nP),C) )
      ENDDO

    ENDIF
    
  ENDDO
  
  RETURN
ENDSUBROUTINE Save_EquPhase

SUBROUTINE Equil_CalcFasAff
!--
!-- calculate affinity of pure phases (i.e. non-aqu'species)
!--
  USE M_Global_Vars,ONLY: vSpc,vFas
  USE M_Basis_Vars, ONLY: tNuFas,vOrdPr
  USE M_Equil_Vars, ONLY: vFasAff
  !
  INTEGER  :: iFs,nFs,nCp
  
  nFs= SIZE(vFas)
  nCp= SIZE(vOrdPr)
  
  IF (nFs>0) THEN
    !note:
    !vDeltaG_Ms(iMs)=  vSpc(vOrdMs(iMs))%G0rt 
    !               - dot_product(tNuMs(iMs,1:nCp),Spc(vOrdPr(1:nCp))%G0rt) 
    !
    !-> deltaG of PRECIPITATION reaction, Prim'Species -> Mineral
    !
    !affinity/RT of FORMATION reaction of mineral,
    != A/RT= -DotProd(Nu,Mu)= -deltaG of Prim'Sp->Mineral, A/RT=ln(K/Q)
    !
    !LnQsK(iMs)= dot_product( tNuMs(iMs,1:nCp), vSpc(vOrdPr(1:nCp))%LnAct ) 
    !          - vDeltaG_Ms(iMs)
    !
    DO iFs=1,nFs
      vFasAff(iFs)= &
      & vFas(iFs)%Grt &
      - DOT_PRODUCT( tNuFas(iFs,1:nCp), &
      &              vSpc(vOrdPr(1:nCp))%Dat%LAct+vSpc(vOrdPr(1:nCp))%G0rt )
      !vFasAff(iFs)= vFasAff(iFs) /SUM(ABS(tNuFas(iFs,1:nCp)))
    ENDDO !_______________________________END
  ENDIF
  
  !expression used in Equil_1:
  !
  !  vFasAff(iFs)= vDeltaG_Fs(iFs) &
  !  &           - dot_product(tNuFas(iFs,1:nCp),vLnAct(vOrdPr(1:nCp)))
  !
  !with
  !
  !  vDeltaG_Fs(:)= vFas(:)%Grt 
  !               - MATMUL(tNuFas(:,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
  !
  !-> can also calculate vFasAff as
  !
  !  vFasAff(iFs)= &
  !  & vSpc(vFas(iFs)%iSpc)%G0rt &
  !  - DOT_PRODUCT(tNuFas(iFs,1:nCp), &
  !  &             vSpc(vOrdPr(1:nCp))%Dat%LAct+vSpc(vOrdPr(1:nCp))%G0rt)
  
  RETURN
ENDSUBROUTINE Equil_CalcFasAff

SUBROUTINE Equil_Trace_Close
!--
!-- close special trace files for activities, Newton, etc
!--
  USE M_Equil_Vars, ONLY: fTrAct, fActiz, fTrcEq
  USE M_Numeric_Tools,ONLY: fNewtF,fNewtR
  !
  IF(fTrAct>0)  THEN; CLOSE(fTrAct); fTrAct= 0; ENDIF
  IF(fActiZ>0)  THEN; CLOSE(fActiZ); fActiZ= 0; ENDIF
  IF(fTrcEq>0)  THEN; CLOSE(fTrcEq); fTrcEq= 0; ENDIF
  !
  IF(fNewtF>0)  THEN; CLOSE(fNewtF); fNewtF= 0; ENDIF
  IF(fNewtR>0)  THEN; CLOSE(fNewtR); fNewtR= 0; ENDIF
  !
  RETURN
ENDSUBROUTINE Equil_Trace_Close

SUBROUTINE Equil_Trace_Init
!--
!-- open special trace files for activities, Newton, etc
!--
  USE M_IOTools,    ONLY: GetUnit
  USE M_Files,      ONLY: DirOut,Files_Index_Write
  USE M_Numeric_Tools,ONLY: fNewtF,fNewtR,fNewt_I
  !
  USE M_Global_Vars,ONLY: vSpc
  USE M_Equil_Write,ONLY: Equil_Write_EnTete
  USE M_Basis_Vars, ONLY: vLCi,vLAs,vLAx,vPrmFw
  USE M_Equil_Vars, ONLY: DebActiv,fTrAct,fActiZ,DebNewt
  !
  INTEGER:: I
  
  IF(DebActiv) THEN
    !
    CALL GetUnit(fTrAct)
    OPEN(fTrAct,FILE=TRIM(DirOut)//"_activ_deb.log")
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirOut)//"_activ_deb.log",&
    & "EQUIL/LOG: activ' aqu'species")
    WRITE(fTrAct,'(A,A1)',ADVANCE='NO') "indx",T_
    DO I=1,SIZE(vSpc)
      IF(vSpc(I)%Typ=="AQU") WRITE(fTrAct,'(A,A1)',ADVANCE='NO') TRIM(vSpc(I)%NamSp),T_
    ENDDO
    WRITE(fTrAct,*)
    !
    CALL GetUnit(fActiz); OPEN(fActiz,FILE=TRIM(DirOut)//"_activ.log")
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirOut)//"_activ.log",&
    & "EQUIL/LOG: log10(Activ)")
    CALL Equil_Write_EnTete(fActiz,vSpc,vPrmFw,vLCi.OR.vLAs.OR.vLAx,Str1="indx")
    !
  ELSE
    !
    fTrAct= 0
    fActiZ= 0
    !
  ENDIF
  
  IF(DebNewt) THEN
    !
    CALL GetUnit(fNewtF)
    OPEN(fNewtF,FILE=TRIM(DirOut)//"_newtequspc1.log")
    CALL Equil_Write_EnTete(fNewtF,vSpc,vPrmFw,vLCi.OR.vLAs,Str1="count1",Str2="count2")
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirOut)//"_newtequspc1.log",&
    & "EQUIL/LOG: Newton solution")
    !
    CALL GetUnit(fNewtR)
    OPEN(fNewtR,FILE=TRIM(DirOut)//"_newtequspc2.log")
    CALL Equil_Write_EnTete(fNewtR,vSpc,vPrmFw,vLCi.OR.vLAs,Str1="count1",Str2="count2")
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirOut)//"_newtequspc2.log",&
    & "EQUIL/LOG: Newton residue")
    !
    fNewt_I=0
    !
  ELSE
    !
    fNewtF=0
    fNewtR=0
    !
  ENDIF
  
  RETURN
ENDSUBROUTINE Equil_Trace_Init

SUBROUTINE Equil_FasTrace_EnTete(vFas,vYesList,F)

  USE M_IOTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut
  USE M_T_Phase,ONLY: T_Phase
  !
  TYPE(T_Phase),INTENT(IN) :: vFas(:)
  LOGICAL,      INTENT(IN) :: vYesList(:)
  INTEGER,      INTENT(OUT):: F
  !
  INTEGER:: iFs,nFs
  
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_equil.log")
  !F=51
  
  WRITE(F,'(3(A,A1))',ADVANCE="NO") "iCount",T_,"nIts",T_,"dXi",T_
  
  nFs=SIZE(vFas)
  !
  DO iFs=1,nFs
    IF(vYesList(iFs)) &
    & WRITE(F,'(A,A1)',ADVANCE="NO") "Mol_"//TRIM(vFas(iFs)%NamFs),T_
  ENDDO
  DO iFs=1,nFs
    IF(vYesList(iFs)) &
    & WRITE(F,'(A,A1)',ADVANCE="NO") "lQsK_"//TRIM(vFas(iFs)%NamFs),T_
  ENDDO
  WRITE(F,*)
  
  RETURN
ENDSUBROUTINE Equil_FasTrace_EnTete

SUBROUTINE Equil_FasTrace_Write(F,ICount,nIts,dXi,vYesList,vFasMole,vFasAff)
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_IOTools,ONLY: OutStrVec
  !
  INTEGER, INTENT(IN):: F
  INTEGER, INTENT(IN):: iCount
  INTEGER, INTENT(IN):: nIts
  REAL(dp),INTENT(IN):: dXi
  LOGICAL, DIMENSION(:),INTENT(IN):: vYesList
  REAL(dp),DIMENSION(:),INTENT(IN):: vFasMole,vFasAff
  
  WRITE(F,'(2(I7,A1),G15.6,A1)',ADVANCE="NO") iCount,T_,nIts,T_,dXi,T_
  CALL OutStrVec(F,vFasMole,     vL=vYesList,CR=.FALSE.)
  CALL OutStrVec(F,-vFasAff/Ln10,vL=vYesList)
  
  RETURN
ENDSUBROUTINE Equil_FasTrace_Write

ENDMODULE M_Equil_Tools

