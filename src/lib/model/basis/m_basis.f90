MODULE M_Basis
!--
!-- basic tools used by higher level modules:
!-- build vFas, vMixModel, ...
!-- operations on basis
!--
  USE M_Kinds
  USE M_Trace,ONLY: fTrc,T_,Stop_,iDebug
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Basis_Build
  PUBLIC:: Basis_Change
  PUBLIC:: Basis_Change_Wrk
  PUBLIC:: Basis_InitBalance
  PUBLIC:: Basis_CpnInert
  !
CONTAINS

SUBROUTINE Basis_Build(LBuffer,vCpn)
  USE M_T_Component,ONLY: T_Component
  USE M_Basis_Tools

  USE M_Global_Vars,ONLY: vEle,vSpc,vFas,vMixFas,vMixModel
  USE M_Global_Vars,ONLY: tFormula
  USE M_Global_Vars,ONLY: nAq,nMn,nGs
  USE M_System_Vars,ONLY: BufferIsExtern,CpnIsSpc

  USE M_Basis_Vars,ONLY: isW,iO_,iH_,iOx,isH_,isOH,isO2,MWsv,iBal
  USE M_Basis_Vars,ONLY: tStoikio
  USE M_Basis_Vars,ONLY: vLCi,vLAx,vLMx,vLAs,vLMs
  USE M_Basis_Vars,ONLY: nCi, nAx, nMx, nAs, nMs
  USE M_Basis_Vars,ONLY: vOrdPr,vOrdAq,vOrdAs,vOrdMs
  USE M_Basis_Vars,ONLY: vPrmFw,vPrmBk
  USE M_Basis_Vars,ONLY: tAlfSp,tAlfPr,tAlfAs,tAlfMs
  USE M_Basis_Vars,ONLY: tNuSp,tNuAs,tNuMs
  USE M_Basis_Vars,ONLY: tNuFas,tAlfFs
  !
  LOGICAL,          INTENT(IN)   :: LBuffer
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  INTEGER:: nSp,nCp,nFs,I
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Build"
  !
  nSp= SIZE(vSpc)
  nCp= SIZE(vEle)
  nFs= SIZE(vFas)
  !
  CALL Basis_Dimensions(   & !
  & BufferIsExtern,        & !IN
  & vCpn,vSpc,nAq,nMn,nGs, & !IN
  & nCi,nMx,nAx,nMs,nAs)     !OUT, may change when components'status change
  !
  CALL Basis_Alloc(nSp,nCp,nFs,nAq,nAs,nMs)
  !
  !---------------------------------------- build logical tables vLxx --
  CALL Basis_SpcTypes_Build( &
  & vSpc,vCpn, &                 !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs)    !OUT
  !
  !----------------------------------- order components, inert/mobile --
  CALL Basis_Cpn_Order( &
  & vSpc,COUNT(vLCi),COUNT(vLAx), & !IN
  & vCpn) !INOUT
  !
  !----------------------------- indexes of elements O, H, Ox in vCpn --
  !--- (may change in basis change, because element'order may change) --
  CALL Basis_Calc_iO_iH_iOx( &
  & vEle,vCpn, &
  & iO_,iH_,iOx)
  !
  !--------- component status and order may have changed -> find iBal --
  IF(ALL(vCpn(:)%Statut=="INERT")) THEN
    iBal=0 !-> no balance element -> mass balance on H
  ELSE
    CALL Basis_FindIBal(vCpn,vSpc,iH_,iBal)
    IF(iBal==0) iBal= iH_
    ! when iBal=iH_, mass balance on H is "replaced" by charge balance
  ENDIF
  !
  !------------ component order may have changed -> new stoikio table --
  tStoikio(1:nCp,1:nSp)= tFormula(vCpn(1:nCp)%iEle,1:nSp) !-> line permutation
  !
  !~ !--- NEW2010-10-29
  !~ IF(CpnIsSpc) THEN
    !~ DO i=1,nCp
      !~ !vCpn(i)%iEle= I
      !~ vCpn(i)%vStoikCp(0:nCp+1)= vSpc(vCpn(i)%iSpc)%vStoikio(0:nCp+1)
    !~ ENDDO
  !~ ELSE
    !~ DO i=1,nCp
      !~ vCpn(i)%iEle= I
      !~ vCpn(i)%vStoikCp(0)= 1 !formula divider !!
      !~ vCpn(i)%vStoikCp(:)= 0
      !~ vCpn(i)%vStoikCp(i)= 1
    !~ ENDDO
  !~ ENDIF
  !~ !---/ NEW2010-10-29
  !
  !---------------------------------------------- build vOrdXx,vPrmXx --
  CALL Basis_SpcOrdPrm_Build(    & !
  & vSpc,vCpn,                   & !IN
  & vLCi,vLAx,vLMx,vLAs,vLMs,    & !IN
  & vOrdPr,vOrdAq,vOrdAs,vOrdMs, & !OUT
  & vPrmFw,vPrmBk)                 !OUT
  !
  IF(iDebug>1) &
  & CALL Basis_StoikTable_Sho(vEle,vCpn,vSpc,tStoikio,vLCi,vLAx,vLMx,vLAs,vLMs)
  !
  !----------------------------------------------- build tAlfXx,tNuXx --
  CALL Basis_AlfaNu_Build( &
  & vEle,vCpn,vSpc,tStoikio,     & !IN
  & vFas,vMixFas,vMixModel,      & !
  & vLCi,vLAx,vLMx,vLAs,vLMs,    & !
  & iDebug>1,LBuffer,            & !
  & vOrdPr,vOrdAs,vOrdMs,        & !
  & vCpn(:)%Mole,                & !INOUT
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !OUT
  & tNuSp,tNuAs,tNuMs,           & !
  & tNuFas,tAlfFs)                 !
  !
  !---------------------------------- modify stoichio tables for iBal --
  IF(iBal/=0) CALL Basis_InitBalance( & !
  & iBal,                 & !IN
  & vOrdPr,vOrdAs,vOrdMs, & !IN
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !INOUT
  & tAlfFs)                 !INOUT
  !
  DO I=1,nCp
    IF(vCpn(I)%Statut/="INERT") vCpn(iBal)%Mole= Zero
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Build"
  !
ENDSUBROUTINE Basis_Build

SUBROUTINE Basis_Change( & !
& Cod,                   & !
!! & tStoikio,           & !
!! & isW,isH_,isOH,isO2, & !
!! & MWsv,               & !
& vCpn)
  
  USE M_T_Component,ONLY: T_Component
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_Basis_Tools
  !
  USE M_Global_Vars,ONLY: vSpc,nAq,vMixFas,vEle !,vSpcAq
  USE M_Basis_Vars, ONLY: isW,isH_,isOH,isO2,tStoikio,MWsv
  !
  CHARACTER(LEN=3), INTENT(IN)   :: Cod
  !! REAL(dp),         INTENT(IN)   :: tStoikio(:,:)
  !! INTEGER,          INTENT(OUT)  :: isW,isH_,isOH,isO2
  !! REAL(dp),         INTENT(OUT)  :: MWsv
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  INTEGER:: nCp,I
  LOGICAL:: LBuffer
  INTEGER,ALLOCATABLE:: vIndex(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Change"
  !
  CALL Basis_IndexInit(isW,isH_,isOH,isO2,MWsv)
  !
  LBuffer= (Cod(1:3)=="DYN") !! .TRUE. !! !.OR. (Cod(1:2)=="EQ")
  !
  nCp= SIZE(vCpn)
  !
  IF(SIZE(vMixFas)>0) &
  & CALL Basis_MixPhase_CpnUpdate( &
  & vMixFas,   & !in
  & vSpc,vCpn)   !inout
  !
  SELECT CASE(TRIM(Cod))
  !
  CASE("DYN","EQU")
    !
    WHERE(vCpn(:)%Statut/="BUFFER") vCpn(:)%Statut="INERT"
    !
    !~ print *,"SUM(tStoikio(I,:) *vSpc(:)%Dat%Mole)"
    !~ DO I=1,nCp
      !~ vCpn(I)%Mole= SUM(tStoikio(I,:) *vSpc(:)%Dat%Mole)
      !~ print '(I3,G15.6,A)',I,vCpn(I)%Mole,TRIM(vEle(vCpn(I)%iEle)%NamEl)
    !~ ENDDO
    !
    ! compute element mole numbers refer to abundances in fluid ONLY
    ! (-> normal case when coming from a speciation calc')
    !print *,"SUM(tStoikio(I,:) *vSpc(:)%Dat%Mole)"
    DO I=1,SIZE(vCpn)
      vCpn(I)%Mole= &
      & SUM(tStoikio(I,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> abundance in fluid
      !print '(I3,G15.6,A)',I,vCpn(I)%Mole,TRIM(vEle(vCpn(I)%iEle)%NamEl)
      !~ & + SUM(tAlfFs(I,1:nFs) *vFas(1:nFs)%Mole)
      !~ WRITE(FF,'(A,A1,G15.6)') vEle(vCpn(I)%iEle)%NamEl,T_,X
    ENDDO
    !~ call pause_
    !
    vCpn(1:nCp)%LnAct= vSpc(vCpn(1:nCp)%iSpc)%Dat%LAct
    !
    !-------------------- assign prim'species to each inert component --
    ALLOCATE(vIndex(nCp))
    !
    CALL Basis_Cpn_ISpc_Select( & !
    ! CALL Basis_Species_Select( & !
    & vSpc,isW,vCpn, & !in
    & vIndex) !inout, only values of vCpn(:)%iSpc are possibly changed
    !
    vCpn(:)%iSpc= vIndex(:)
    DEALLOCATE(vIndex)
    !-------------------/ assign prim'species to each inert component --
    !
    !! IF(iOx/=0) CALL Basis_RedoxAssign(vCpn,vEle,vSpc)
  END SELECT
  !
  CALL Basis_Build(LBuffer,vCpn)
  !-> include following
  ! Basis_Dimensions      !-> nCi,nMx,nAx,nMs,nAs
  ! Basis_Alloc           !-> allocate basis_vars
  ! Basis_SpcTypes_Build  !-> build vLCi,vLAx,vLMx,vLAs,vLMs
  ! Basis_Cpn_Order       !-> order components, inert/mobile

  ! Basis_Calc_iO_iH_iOx  !-> indexes of elements O, H, Ox in vCpn
  ! Basis_SpcOrdPrm_Build !-> vOrdPr,vOrdAq,vOrdAs,vOrdMs,vPrmFw,vPrmBk
  ! Basis_AlfaNu_Build    !-> build tAlfXx,tNuXx
  ! Basis_InitBalance     !
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Change"
  !
END SUBROUTINE Basis_Change

SUBROUTINE Basis_Change_Wrk(isW,vCpn)
  !---------------------------------------------------------------------
  USE M_T_Component,ONLY: T_Component
  USE M_Basis_Tools,ONLY: Basis_Cpn_ISpc_Select
  !
  USE M_Global_Vars,ONLY: vSpc
  !---------------------------------------------------------------------
  INTEGER,          INTENT(IN)   :: isW
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  INTEGER,ALLOCATABLE:: vIndex(:)
  !
  IF(iDebug==4) WRITE(fTrc,'(/,A)') "< Basis_Change_wrk"
  !
  ALLOCATE(vIndex(SIZE(vCpn)))
  !
  CALL Basis_Cpn_ISpc_Select( & !
  ! CALL Basis_Species_Select( & !
  & vSpc,isW,vCpn, & !in
  & vIndex)
  !
  !--------- change basis only if the primary species set has changed --
  IF(ANY (vCpn(:)%iSpc /= vIndex(:)) ) THEN
    !
    vCpn(:)%iSpc= vIndex(:)
    !
    CALL Basis_Build(.FALSE.,vCpn)
    !
      IF(iDebug>2) WRITE(*,'(A)') "===============< basis changed >===="
    !
  ENDIF
  !
  DEALLOCATE(vIndex)
  !
  IF(iDebug==4) WRITE(fTrc,'(A,/)') "</ Basis_Change_wrk"
END SUBROUTINE Basis_Change_Wrk

SUBROUTINE Basis_IndexInit(isW,isH_,isOH,isO2,MWsv)
!--
!-- calc' dimensions and indexes that are constant for given vSpc
!--
  USE M_T_Species,  ONLY: Species_Index
  !
  USE M_Global_Vars,ONLY: vSpc,nAq,nMn,nGs
  !
  INTEGER, INTENT(OUT):: isW,isH_,isOH,isO2
  REAL(dp),INTENT(OUT):: MWsv
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_IndexInit"
  !
  nAq= COUNT(vSpc%Typ(1:3)=="AQU")
  nMn= COUNT(vSpc%Typ(1:3)=="MIN")
  nGs= COUNT(vSpc%Typ(1:3)=="GAS")
  nMn= nMn + nGs
  !
  ! find indexes of solvent, H+, OH-, O2 in vSpc
  ! -> unchanged in basis switch -> called only once for a given vSpc
  !-------------- ranks of species H+ and OH- in current species list --
  isW=  Species_Index("H2O",  vSpc)
  isOH= Species_Index("OH-",  vSpc); IF(isOH==0) isOH=Species_Index("OH[-]",vSpc)
  isH_= Species_Index("H+",   vSpc); IF(isH_==0) isH_=Species_Index("H[+]",vSpc)
  isO2= Species_Index("O2,AQ",vSpc); IF(isO2==0) isO2=Species_Index("O2(AQ)",vSpc)
  isO2= Species_Index("O2_AQ",vSpc)
  !
  IF(isW==0)  CALL Stop_("species H2O not found !!!") !=============STOP
  IF(isH_==0) CALL Stop_("species H+ not found  !!!") !=============STOP
  IF(isOH==0) CALL Stop_("species OH- not found !!!") !=============STOP
  !
  MWsv=vSpc(isW)%WeitKg
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_IndexInit"
  !
ENDSUBROUTINE Basis_IndexInit

SUBROUTINE Basis_InitBalance( & !
& iBal,                 & !IN
& vOrdPr,vOrdAs,vOrdMs, & !IN
& tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs) !INOUT
!--
!-- replace current stoichio value by the species charge
!-- in line iBal of all stoikio matrices tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs
!--
  !
  USE M_Global_Vars,ONLY: vSpc
  !
  INTEGER, INTENT(IN)   :: iBal
  INTEGER, INTENT(IN)   :: vOrdPr(:),vOrdAs(:),vOrdMs(:)
  REAL(dp),INTENT(INOUT):: tAlfSp(:,:),tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:),tAlfFs(:,:)
  !
  tAlfSp(iBal,:)=                    REAL(vSpc(:)%Z)
  tAlfPr(iBal,:)=                    REAL(vSpc(vOrdPr(:))%Z)
  IF(SIZE(vOrdAs)>0) tAlfAs(iBal,:)= REAL(vSpc(vOrdAs(:))%Z)
  IF(SIZE(vOrdMs)>0) tAlfMs(iBal,:)= REAL(vSpc(vOrdMs(:))%Z)
  tAlfFs(iBal,:)= Zero
  !
ENDSUBROUTINE Basis_InitBalance

SUBROUTINE Basis_CpnInert(isW,tStoikio,vCpnIn,vCpnInert)
!--
!-- used only in Equil_Write,
!-- for printing composition of aqueous phase as closed system
!-- make "INERT" all "MOBILE" components,
!-- reorder components, assign prim'species, ...
!--
  USE M_T_Component,ONLY: T_Component
  USE M_Basis_Tools,ONLY: Basis_Species_Select_Test
  USE M_Basis_Tools,ONLY: Basis_Cpn_ISpc_Select
  !
  USE M_Global_Vars,ONLY: vEle,vSpc,nAq,tFormula
  !
  INTEGER,          INTENT(IN) :: isW
  REAL(dp),         INTENT(IN) :: tStoikio(:,:)
  TYPE(T_Component),INTENT(IN) :: vCpnIn(:)
  TYPE(T_Component),INTENT(OUT):: vCpnInert(:)
  !
  INTEGER:: i
  INTEGER,ALLOCATABLE:: vIndex(:)
  !
  vCpnInert= vCpnIn
  vCpnInert(:)%Statut="INERT"
  !
  !-------- compute total amount of each component in the fluid phase --
  DO i=1,SIZE(vCpnInert)
    vCpnInert(i)%Mole= &
    & SUM( tStoikio(i,:)*vSpc(:)%Dat%Mole, MASK= (vSpc(:)%Typ=="AQU") )
  ENDDO
  !
  !---------------------------------------------- assign prim'species --
  ALLOCATE(vIndex(SIZE(vCpnInert)))
  !
  !~ CALL Basis_Species_Select_Test( & !
  !~ & vSpc,isW,vCpnInert,tFormula, & !in
  !~ & vIndex)

  CALL Basis_Cpn_ISpc_Select( & !
  ! CALL Basis_Species_Select( & !
  & vSpc,isW,vCpnInert, & !in
  & vIndex) !out
  !
  vCpnInert(:)%iSpc= vIndex(:)
  !---------------------------------------------/ assign prim'species --
  !
  DEALLOCATE(vIndex)
  !
ENDSUBROUTINE Basis_CpnInert

ENDMODULE M_Basis


