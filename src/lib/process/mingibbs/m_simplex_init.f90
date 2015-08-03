MODULE M_Simplex_Build
!--
!-- io module for simplex calculation
!-- independant (provisionally) from other io modules
!--
  USE M_Kinds
  USE M_Trace,      ONLY: fTrc,iDebug,T_,Stop_,Pause_,Warning_
  USE M_T_Component,ONLY: T_Component, Component_Zero
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Simplex_Build
  PUBLIC:: Simplex_CloseFiles
  !
  INTEGER:: nFasPur,nFasSol
  INTEGER:: Discret
  != dimension of discretization (=number of points on each parameter)
  !
  !! LOGICAL:: Redox
  !
  INTEGER,PUBLIC:: fSpl1=0, fSpl2=0
  ! files that should remain accessible
  ! by different subroutines throughout program
  !
CONTAINS

SUBROUTINE Simplex_Build
  !
  USE M_T_Element,   ONLY: Element_Index
  USE M_T_Species,   ONLY: T_Species,Species_Stoikio_Calc
  USE M_Dtb_Const,   ONLY: Tref,Pref,T_CK
  !
  USE M_Global_Alloc,ONLY: MixModels_Alloc,Phases_Alloc,MixPhases_Alloc
  USE M_Global_Alloc,ONLY: DiscretParam_Alloc
  USE M_Global_Tools,ONLY: Global_TP_Update
  !
  USE M_DiscretModel_Read
  USE M_DiscretModel_Tools !DiscretParam_Init
  !
  USE M_Global_Vars, ONLY: vEle,vSpcDtb,vSpc,vMixModel,vMixFas,vFas
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam
  !
  USE M_GEM_Vars, ONLY: TdgK,Pbar,vCpnGEM,tStoikioGEM
  !---------------------------------------------------------------------
  LOGICAL:: Ok
  INTEGER:: N
  TYPE(T_Species),ALLOCATABLE:: vSpcTmp(:)
  !---------------------------------------------------------------------
  IF(iDebug>2) WRITE(fTrc,'(/,A)') "<---------------------Simplex_Build"
  !
  nFasPur=0
  nFasSol=0
  !
  TdgK= Tref
  Pbar= Pref
  !
  ! IF(Element_Index("O__",vEle)==0) &
  ! & CALL Stop_("Element Oxygen not found in vEle")
  !
  !--build vCpnGEM
  CALL Component_Read(vEle)
  !
  !--rebuild species array vSpc consistent with vCpnGEM
  CALL Species_Select(SIZE(vEle),vCpnGEM)
  !
  !--build new vEle
  CALL Element_Select(vCpnGEM)
  !
  CALL Element_Sort(vEle)
  !
  CALL Component_Stoikio_Calc(vEle,vCpnGEM,Ok)
  !
  !--build species array vSpc consistent with vEle
  !--OBSOLETE
  ! CALL Spl_Species_Alloc(vEle)
  !
  CALL Species_Stoikio_Calc(vEle,vSpc,Ok)
  !
  !--build MixModel database -> vMixModel
  CALL MixModels_Alloc(vSpc)
  !
  !-------------------------------------------append discretized species
  !--read discretization models, build vDiscretModel
  CALL DiscretModel_Read(vMixModel)
  !
  N= SIZE(vDiscretModel)
  !
  IF(N>0) THEN
    !
    CALL DiscretParam_Alloc(vDiscretModel) !-> allocate vDiscretParam
    !
    ALLOCATE(vSpcTmp(SIZE(vDiscretParam)))
    !
    CALL DiscretParam_Init( &
    & vEle,vSpc,vMixModel,vDiscretModel, &
    & vDiscretParam,vSpcTmp) !-> build vSpcDiscret
    !
    CALL Species_Append(vSpcTmp) !-> new vSpc !!!
    !
    DEALLOCATE(vSpcTmp)
    !
    CALL DiscretSpecies_Stoikio_Calc( & !
    & vEle,          & !IN
    & vMixModel,     & !IN
    & vDiscretModel, & !IN
    & vDiscretParam, & !IN
    & vSpc)            !INOUT
    !
  ENDIF
  !------------------------------------------/append discretized species
  !
  !--- read phase compositions, build vMixFas
  CALL MixPhases_Alloc(vSpc,vMixModel)
  !
  !--- build vFas
  CALL Phases_Alloc(vSpc,vMixFas)
  !
  IF(iDebug>2) PRINT *,"TdgC/Pbar=", TdgK-T_CK,Pbar  !;  pause
  !
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  !
  !--- compute tFormula
  CALL Spl_Tables_Build
  !
  !--- allocate, compute tStoikioGEM
  IF(ALLOCATED(tStoikioGEM)) DEALLOCATE(tStoikioGEM)
  ALLOCATE(tStoikioGEM(SIZE(vFas),SIZE(vCpnGEM)))
  !
  CALL Spl_Stoikio_Calc(vCpnGEM,tStoikioGEM)
  !---/
  !
  !--- write system stoikiometry
  IF(iDebug>2) CALL Spl_WriteSystem(vCpnGEM,tStoikioGEM)
  !
  IF(iDebug>2) WRITE(fTrc,'(A,/)') "</--------------------Simplex_Build"
  !
  RETURN
ENDSUBROUTINE Simplex_Build

SUBROUTINE Species_Select(nEl,vCpn)
!--
!-- rebuild species array vSpc consistent with vCpnGEM
!--
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  !
  USE M_Global_Vars, ONLY: vSpc
  !
  INTEGER, INTENT(IN):: nEl
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  TYPE(T_Species),ALLOCATABLE:: vTmp(:)
  REAL(dp),ALLOCATABLE:: tStoikio(:,:) !tStoikio(:,:),
  INTEGER :: I,iS,nS,nC,N
  LOGICAL :: Ok

  IF(iDebug>2) WRITE(fTrc,'(/,A)') "< Species_Select"
  
  nC= SIZE(vCpn)
  nS= SIZE(vSpc)
  ALLOCATE(vTmp(nS))
  vTmp(:)= vSpc(:)

  ALLOCATE(tStoikio(nC+1,nEl))
  DO I=1,nC
    tStoikio(I,1:nEl)= vCpn(I)%vStoikCp(1:nEl)  !tStoikio(I,:)
  ENDDO

  n= 0
  DO iS=1,nS
  
    IF(vSpc(iS)%Typ /= "AQU" .OR. TRIM(vSpc(iS)%NamSp)=="H2O") THEN
    
      tStoikio(nC+1,1:nEl)= vSpc(iS)%vStoikio(1:nEl)
      
      CALL Check_Independent(tStoikio,Ok)
      
      IF(.NOT. Ok) THEN
        ! if species is dependent on basis vCpn(:),
        ! then add to vSpc
        N= N+1
        vTmp(N)= vSpc(iS)
        !~ WRITE(6,'(2A)') "OK SPECIES= ", TRIM(vSpc(is)%NamSp)
        IF(iDebug>2) WRITE(fTrc,'(A)') TRIM(vSpc(iS)%NamSp)
      ENDIF
      
    ENDIF
  
  ENDDO

  DEALLOCATE(vSpc)
  ALLOCATE(vSpc(N))
  vSpc(1:N)= vTmp(1:N)

  DEALLOCATE(vTmp)
  DEALLOCATE(tStoikio)

  IF(iDebug>2) WRITE(fTrc,'(A,/)') "</ Species_Select"
  
  RETURN
END SUBROUTINE Species_Select

SUBROUTINE Check_Independent(tStoikio,Ok)
  USE M_Numeric_Tools,ONLY: iMinLoc_R,iFirstLoc,iMaxLoc_R
  !
  REAL(dp),INTENT(IN) :: tStoikio(:,:)
  LOGICAL, INTENT(OUT):: Ok
  !
  REAL(dp),ALLOCATABLE:: tStoikTmp(:,:) !tStoikio(:,:),
  INTEGER, ALLOCATABLE:: vIPivot(:)
  !
  INTEGER :: I,J,K
  INTEGER :: N,nEl
  REAL(dp):: Y,Pivot
  
  N=   SIZE(tStoikio,1)
  nEl= SIZE(tStoikio,2)
  ALLOCATE(tStoikTmp(N,nEl))     ;  tStoikTmp(:,:)= Zero
  ALLOCATE(vIPivot(nEl))         ;  vIPivot(:)= 0
  !
  Ok= .TRUE.
  K= 0
  DO I=1,N
    !
    K= K+1
    tStoikTmp(K,:)= tStoikio(I,:)
    DO J=1,K-1
      Y= tStoikTmp(K,vIPivot(J))
      tStoikTmp(K,:)= tStoikTmp(K,:)- Y*tStoikTmp(J,:)
    ENDDO
    vIPivot(K)= iMaxLoc_R(ABS(tStoikTmp(K,:)))
    Pivot= tStoikTmp(K,vIPivot(K))
    IF(ABS(Pivot)<1.D-9) THEN
      != all coeff's are 0 -> this species is not independent -> EXIT ==
      Ok= .FALSE.
      DEALLOCATE(tStoikTmp,vIPivot) !tStoikio,
      RETURN
    ENDIF
    !
    tStoikTmp(K,:)= tStoikTmp(K,:) /Pivot
    !
  ENDDO
  !
  ! WRITE(sFMT,'(a,i3,a)') '(',nEl,'(G12.3,1X))'
  ! DO I=1,N
  !   WRITE(12,sFMT) (tStoikTmp(I,K),K=1,nEl)
  ! ENDDO
  !
  DEALLOCATE(tStoikTmp,vIPivot) !tStoikio,
  !
END SUBROUTINE Check_Independent

SUBROUTINE Species_Append(vSpcAdd)
!--
!-- append vSpcAdd to current vSpc
!-- -> produce a new vSpc
!--
  USE M_T_Species, ONLY: T_Species
  !
  USE M_Global_Vars, ONLY: vSpc
  !
  TYPE(T_Species),INTENT(IN):: vSpcAdd(:)
  !
  TYPE(T_Species),ALLOCATABLE:: vSpcAll(:)
  INTEGER:: M, N
  !
  M= SIZE(vSpc)
  N= SIZE(vSpcAdd)
  !
  ALLOCATE(vSpcAll(M+N))
  vSpcAll(1   :M )= vSpc(1:M)
  vSpcAll(M+1:M+N)= vSpcAdd(1:N)
  !
  DEALLOCATE(vSpc)
  ALLOCATE(vSpc(N+M)) ; vSpc= vSpcAll
  !
  DEALLOCATE(vSpcAll)
  !
ENDSUBROUTINE Species_Append
  !
SUBROUTINE Simplex_CloseFiles
  IF(fSpl1>0) CLOSE(fSpl1) ; fSpl1= 0
  IF(fSpl2>0) CLOSE(fSpl2) ; fSpl2= 0
ENDSUBROUTINE Simplex_CloseFiles

SUBROUTINE Component_Read(vEle)
!--
!-- read components from SYSTEM.GEM block -> build vCpnGEM
!--
  USE M_IOTools,  ONLY: Str_Append
  USE M_Files,    ONLY: NamFInn
  USE M_Dtb_Const,ONLY: T_CK
  USE M_IoTools,  ONLY: GetUnit,LinToWrd,AppendToEnd,WrdToReal,Str_Upper
  USE M_T_Element,ONLY: T_Element,Formula_Read,Element_Index,Formula_Build
  USE M_Dtb_Read_Tools
  !
  USE M_GEM_Vars,ONLY: vCpnGEM,TdgK,Pbar
  !
  TYPE(T_Element),INTENT(IN):: vEle(:)
  !
  CHARACTER(LEN=255)    :: L,W
  TYPE(T_Component)     :: Cpn
  LOGICAL :: sOk,fOk,CpnOk
  INTEGER :: nEl,N,i,iOx
  INTEGER :: ZSp,nDiv,Ztest
  INTEGER :: f,ios
  LOGICAL :: UnitIsMole,UnitIsGram
  REAL(dp):: X
  CHARACTER(LEN=4):: Str
  !
  REAL(dp),ALLOCATABLE:: tStoikio(:,:)
  INTEGER, ALLOCATABLE:: vStoik(:)
  TYPE(T_Component),ALLOCATABLE:: vC(:)
  !
  LOGICAL:: EcformIsOk
  INTEGER:: fFormula
  CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE:: vElement
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Component_Read"
  !
  !CodFormula= "ECFORM"
  fFormula= 0
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  iOx= Element_Index("OX_",vEle)
  !
  nEl=SIZE(vEle)
  ALLOCATE(vStoik(1:nEl))
  ALLOCATE(vC(1:nEl))
  !
  IF(iDebug>2) PRINT *,"Component_Build"
  !
  N=0
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,sOk)
    IF(W(1:1)=='!') CYCLE DoFile
    CALL AppendToEnd(L,W,sOk)
    !
    SELECT CASE(W)
    !
    CASE("SYSTEM.SIMPLEX","SYSTEM.GEM")
      !
      DoBlock: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,sOk)
        IF(W(1:1)=='!') CYCLE DoBlock
        CALL AppendToEnd(L,W,sOk)
        !
        UnitIsMole= .TRUE.
        UnitIsGram= .FALSE.
        !
        SELECT CASE(W)
        
        CASE("ENDINPUT")                 ;  EXIT DoFile
        CASE("END","ENDSYSTEM.SIMPLEX")  ;  EXIT DoBlock
        CASE("ENDSYSTEM.GEM")            ;  EXIT DoBlock
        
        CASE("TDGC")
          CALL LinToWrd(L,W,sOk)
          CALL WrdToReal(W,X)
          TdgK= X +T_CK
          CYCLE DoBlock !----------------------------------cycle DoBlock
        
        CASE("TDGK")
          CALL LinToWrd(L,W,sOk)
          CALL WrdToReal(W,X)
          TdgK= X
          CYCLE DoBlock !----------------------------------cycle DoBlock
        
        CASE("PBAR")
          CALL LinToWrd(L,W,sOk)
          CALL WrdToReal(W,X)
          Pbar= X
          CYCLE DoBlock !----------------------------------cycle DoBlock
        
        CASE("MOLE")
          UnitIsMole= .TRUE.
          CALL LinToWrd(L,W,sOk)
        
        CASE("GRAM")
          UnitIsGram= .TRUE.
          CALL LinToWrd(L,W,sOk)
          
        ENDSELECT
        !
        CALL Component_Zero(Cpn)
        !
        CALL Str_Append(W,3)
        Cpn%NamCp= TRIM(W)
        !
        CALL LinToWrd(L,W,sOk,"NO")
        !-> compact formula, character case should be conserved :!!
        !
        !IF(CodFormula=="SCFORM") THEN
          CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
          IF(.NOT.EcformIsOk) THEN
            IF(iDebug>0) PRINT '(3A)',TRIM(W)," =ERROR IN FORMULA, SKIPPED !!"
            CYCLE DoBlock
          ENDIF
        !ENDIF
        !
        CALL Str_Upper(W)
        !
        Cpn%Formula= TRIM(W)
        !
        !----------------------------------------------- process formula 
        !--------------------- if formula is Ok, save component in vC(:) 
        CALL Formula_Read(Cpn%Formula,vEle,ZSp,nDiv,fOk,vStoik)
        !
        Ztest= DOT_PRODUCT(vStoik(1:nEl),vEle(1:nEl)%Z)
        ! IF(fOk .AND. ZSp==Ztest) THEN
        IF(fOk) THEN
          !
          ! vStoikio has globally fixed SIZE, nElMax,
          ! vStoik has local SIZE nEl !!
          !
          Cpn%Statut= "INERT"
          Cpn%vStoikCp(1:nEl)= vStoik(1:nEl)
          Cpn%vStoikCp(0)=     nDiv
          Cpn%vStoikCp(nEl+1)= ZSp
          !
          IF(Zsp/=Ztest) THEN
            IF(iOx/=0) THEN
              Cpn%vStoikCp(iOx)= Zsp - Ztest
              WRITE(Str,'(I4)') Cpn%vStoikCp(iOx)
              Cpn%Formula= TRIM(Cpn%Formula)//"OX("//TRIM(ADJUSTL(Str))//")"
            ELSE
              PRINT *,"ZSp,Stoikio=",ZSp,Ztest
              IF(iDebug>0) WRITE(fTrc,'(A,A1,A15,A1,A39)') &
              & "REJECT",T_,Cpn%NamCp,T_,Cpn%Formula
              IF(iDebug>0) PRINT '(3A)',"rejected: ",Cpn%NamCp,Cpn%Formula
              CYCLE
            END IF
          END IF
          !
          !CALL Formula_Build(vEle,Cpn%vStoikCp,Zsp,nDiv,Cpn%Formula)
          !
          Cpn%iMix= 0
          !
          IF(iDebug>0) WRITE(fTrc,'(A,A1,A15,A1,A39,A1,I3)') &
          & "ACCEPT",T_,Cpn%NamCp,T_,Cpn%Formula,T_,N
          IF(iDebug>2) PRINT '(3A)',"accepted: ",Cpn%NamCp,Cpn%Formula
          !
          CALL LinToWrd(L,W,sOk)
          CALL WrdToReal(W,Cpn%Mole)
          IF(UnitIsGram) THEN
            X= DOT_PRODUCT( &
            & Cpn%vStoikCp(1:nEl), &
            & vEle(1:nEl)%WeitKg)/Cpn%vStoikCp(0)
            Cpn%Factor= X *1.0D3
            Cpn%Mole= Cpn%Mole /Cpn%Factor
            !! print *,"Cpn%Factor",Cpn%Factor
          ENDIF
          !
          N=N+1
          vC(N)= Cpn
          !
        ELSE
          !
          PRINT *,"ZSp,Stoikio=",ZSp,Ztest
          IF(iDebug>0) WRITE(fTrc,'(A,A1,A15,A1,A39)') &
          & "REJECT",T_,Cpn%NamCp,T_,Cpn%Formula
          IF(iDebug>0) PRINT '(3A)',"rejected: ",Cpn%NamCp,Cpn%Formula
          !
        ENDIF
        !-------------------------------------------/ process formula --
        !
      ENDDO DoBlock
      !
      IF(N==0) CALL Stop_("< Found No Components...Stop") !===== stop ==
      !
    ENDSELECT
    !
  ENDDO DoFile
  !
  !! IF(iDebug>2) CALL Pause_
  CLOSE(f)
  !
  IF(ALLOCATED(vCpnGEM)) DEALLOCATE(vCpnGEM)
  ALLOCATE(vCpnGEM(1:N))
  !
  vCpnGEM(1:N)= vC(1:N)
  !
  DEALLOCATE(vStoik)
  DEALLOCATE(vC)
  !
  !----------------------------------------------- check independency --
  ALLOCATE(tStoikio(N,nEl))
  DO I=1,N
    tStoikio(I,1:nEl)= vCpnGEM(I)%vStoikCp(1:nEl)  !tStoikio(I,:)
  ENDDO
  CALL Check_Independent(tStoikio,CpnOk)
  IF(.NOT. CpnOk) CALL Stop_("Components Not Independent")
  DEALLOCATE(tStoikio)
  !-----------------------------------------------/check independency --
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Component_Read"
  !
  RETURN
ENDSUBROUTINE Component_Read

SUBROUTINE Element_Select(vCpn)
!--
!-- reduce element list to those involved in vCpn,
!-- then re-build vEle
!--
  USE M_T_Component,ONLY: T_Component
  USE M_T_Element,  ONLY: T_Element
  !
  USE M_Global_Vars,ONLY: vEle
  !
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  INTEGER,ALLOCATABLE:: tStoik(:,:)
  INTEGER:: nEl, i, N
  TYPE(T_Element),ALLOCATABLE:: vE(:)
  !
  IF(iDebug>2) WRITE(fTrc,'(A)') "<----------------------Element_Select"
  IF(iDebug>2) WRITE(fTrc,'(A)') "-> New Element List"
  !
  nEl= SIZE(vEle)
  !
  ALLOCATE(tStoik(SIZE(vCpn),nEl))
  DO i=1,SIZE(vCpn)
    tStoik(i,1:nEl)= vCpn(i)%vStoikCp(1:nEl)
  ENDDO
  !
  IF(iDebug>2) THEN
    WRITE(fTrc,'(*(A,A1))') (vEle(i)%NamEl,t_,i=1,nEl)
    CALL WriteMat_I(fTrc,tStoik)
  END IF
  !
  ALLOCATE(vE(1:nEl))
  N=0
  !
  !---------------------------select elements involved in the components
  DO I=1,nEl
    IF(ANY(tStoik(:,i)/=0)) THEN
      N=N+1
      vE(N)= vEle(I) !-> new element list
      IF(iDebug>2) WRITE(fTrc,'(I3,1X,A)') N,vEle(I)%NamEl
    ENDIF
  ENDDO
  !---/
  !
  DEALLOCATE(vEle); ALLOCATE(vEle(1:N))
  vEle(1:N)= vE(1:N)
  !
  DEALLOCATE(tStoik)
  DEALLOCATE(vE)
  !
  IF(iDebug>2) WRITE(fTrc,'(A,/)') "</-------------------Element_Select"
  !
  RETURN
END SUBROUTINE Element_Select

SUBROUTINE Element_Sort(vEle)
  USE M_T_Element,ONLY: T_Element,Element_Index
  !
  TYPE(T_Element),INTENT(INOUT):: vEle(:)
  !
  TYPE(T_Element):: E
  INTEGER:: i
  !
  i= Element_Index("O__",vEle)
  IF(i/=1) THEN
    E= vEle(1)
    vEle(1)= vEle(i)
    vEle(i)= E
  ENDIF
  !
ENDSUBROUTINE Element_Sort

SUBROUTINE Component_Stoikio_Calc(vEle,vCpn,Ok)
!--
!-- update vCpn(:)%vStoikCp
!--
  USE M_T_Element,  ONLY: T_Element,Element_Index
  USE M_T_Component,ONLY: T_Component,Component_Stoikio
  !
  TYPE(T_Element),  INTENT(IN)   :: vEle(:)
  TYPE(T_Component),INTENT(INOUT):: vCpn(:) !modIFs in nDiv,Z,Oxy
  LOGICAL,          INTENT(OUT)  :: Ok
  !
  INTEGER:: I,J,ieOx
  LOGICAL:: fOk
  !
  IF(iDebug>2)  WRITE(fTrc,'(A)') "< Component_Stoikio_Calc"
  !
  ieOx= Element_Index("OX_",vEle)
  !
  Ok= .true.
  DO I=1,SIZE(vCpn)
    CALL Component_Stoikio(vEle,ieOx,vCpn(I),fOk)
    IF(.NOT. fOk) THEN
      Ok= .false.
      !! Msg= "Stoikiometry problem in species "//TRIM(vSpc(I)%NamSp)
      RETURN
    ENDIF
  ENDDO
  !
  !------------------------------------------------------------ trace --
  IF(iDebug>2) THEN
    DO i=1,SIZE(vEle)
      WRITE(fTrc,'(A,A1)',ADVANCE="no") vEle(i)%NamEl,t_
    ENDDO
    WRITE(fTrc,*)
    DO i=1,SIZE(vCpn)
      DO j=1,SIZE(vEle)
        WRITE(fTrc,'(I3,A1)',ADVANCE="no") vCpn(i)%vStoikCp(j),t_
      ENDDO
      WRITE(fTrc,*)
    ENDDO
  ENDIF
  !-----------------------------------------------------------/ trace --
  IF(iDebug>2)  WRITE(fTrc,'(A)') "</ Component_Stoikio_Calc"
  !
  RETURN
ENDSUBROUTINE Component_Stoikio_Calc

INTEGER FUNCTION iMaxLoc_I(ARR) !"NR"
  REAL(dp),INTENT(IN) :: ARR(:)
  INTEGER:: IMAX(1)
  IMAX=MAXLOC(ARR(:))
  iMaxLoc_I=IMAX(1)
END FUNCTION iMaxLoc_I

SUBROUTINE Spl_Species_Alloc(vEle)
!--
!-- retrieve all species consistent with current vEle
!-- OBSOLETE: replaced by Species_Select
!--
  USE M_T_Element,   ONLY: T_Element,Element_Index
  USE M_T_Species,   ONLY: T_Species,Species_Stoikio
  USE M_Global_Vars, ONLY: vSpcDtb,vSpc
  USE M_Dtb_Calc,    ONLY: DtbSpc_GrtTable_Build
  !
  TYPE(T_Element),INTENT(IN):: vEle(:)
  !
  TYPE(T_Species),ALLOCATABLE:: vSpc_(:)
  INTEGER:: nEl,nSp,ieOx,I,N
  LOGICAL:: fOk
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< SplSpecies_Alloc"
  !
  nSp=  SIZE(vSpc)
  nEl=  SIZE(vEle)
  ieOx= Element_Index("OX_",vEle)
  !
  ALLOCATE(vSpc_(1:nSp))
  !
  !--- select species consistent with vEle
  N=0
  DO I=1,nSp
    IF(vSpc(I)%Typ /= "AQU" .OR. TRIM(vSpc(I)%NamSp)=="H2O") THEN
      CALL Species_Stoikio(vEle,ieOx,vSpc(I),fOk)
      IF(fOk) THEN
        N=N+1
        vSpc_(N)= vSpc(I)
      ENDIF
    ENDIF
  ENDDO
  !
  DEALLOCATE(vSpc); ALLOCATE(vSpc(1:N))
  vSpc(1:N)=   vSpc_(1:N)
  DEALLOCATE(vSpc_)
  !
  !~ nSp=  SIZE(vSpc)
  !~ CALL Warning_("vTPpath must be ALLOCATED before Spl_Species_Alloc")
  !~ nTP=  SIZE(vTPpath)
  !~ IF(ALLOCATED(tGrt)) DEALLOCATE(tGrt); ALLOCATE(tGrt(1:nSp,1:nTP))
  !~ CALL DtbSpc_GrtTable_Build(vTPpath,vSpcDtb,vSpc,tGrt)
  !~ !
  IF(iDebug>0) THEN
    WRITE(fTrc,'(A)') "-> New Species List"
    DO i=1,SIZE(vSpc)
      WRITE(fTrc,'(I3,1X,A)') I,vSpc(I)%NamSp
    ENDDO
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ SplSpecies_Alloc"
  !
ENDSUBROUTINE Spl_Species_Alloc

SUBROUTINE Spl_Tables_Build
!--
!-- build tFormula(1:nFs,1:nEl): table of species stoikios in lines
!--
  USE M_T_Element,   ONLY: Formula_Read
  USE M_T_Species,   ONLY: T_Species
  USE M_Global_Vars, ONLY: vEle,vFas,vSpc,tFormula
  !
  INTEGER:: nEl,nFs,I
  INTEGER,DIMENSION(:),ALLOCATABLE:: vStoik
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Spl_BuildTables"
  !
  nEl=SIZE(vEle)
  nFs=SIZE(vFas)
  !
  DEALLOCATE(tFormula) ;  ALLOCATE(tFormula(1:nFs,1:nEl))
  !
  ALLOCATE(vStoik(0:nEl))
  !
  DO I=1,nFs
    vStoik(0:nEl)= vSpc(vFas(I)%iSpc)%vStoikio(0:nEl)
    tFormula(I,1:nEl)= vStoik(1:nEl) /REAL(vStoik(0))
  ENDDO
  !
  DEALLOCATE(vStoik)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Spl_BuildTables"
ENDSUBROUTINE Spl_Tables_Build

SUBROUTINE Spl_Stoikio_Calc(vCpn,tStoikio)
!--
!-- Compute tStoikioCpn from tStoikEle
!-- tStoikioCpn(1:nFs,1:nCp)- stoikio of phases in terms of components
!--

  USE M_T_Component,ONLY: T_Component
  USE M_T_Element,   ONLY: Formula_Read
  USE M_Numeric_Mat, ONLY: LU_Decomp, LU_BakSub
  !
  USE M_Global_Vars, ONLY: tFormula,vFas,vEle
  !
  TYPE(T_Component),INTENT(IN) :: vCpn(:)
  REAL(dp),         INTENT(OUT):: tStoikio(:,:)
  !
  REAL(dp),ALLOCATABLE:: A(:,:)
  INTEGER, ALLOCATABLE:: Indx(:)
  REAL(dp),ALLOCATABLE:: Y(:)

  ! stoikiometry table of the components in terms of the elements
  REAL(dp),ALLOCATABLE:: tStoik(:,:)
  !
  LOGICAL :: bSingul
  REAL(dp):: D
  INTEGER :: I,nEl,nCp,nFs
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Spl_Stoikio_Calc"
  !
  nFs=SIZE(vFas)
  nCp=SIZE(vCpn)
  nEl=SIZE(vEle)
  !
  ALLOCATE(tStoik(nCp,nEl))
  DO i=1,nCp
    tStoik(i,1:nEl)= vCpn(i)%vStoikCp(1:nEl) /REAL(vCpn(i)%vStoikCp(0))
  ENDDO
  !
  ALLOCATE(A(1:nCp,1:nCp))
  A(1:nCp,1:nCp)= TRANSPOSE(tStoik(1:nCp,2:nEl))
  ! column_1=Oxygen, -> already taken in account by Spl_ReadFormula
  ! -> columns 2:nEl of ElementFormula Matrix gives ElementStoikio of components 1:nCp
  !
  ALLOCATE(Indx(1:nCp))
  CALL LU_Decomp(A,Indx,D,bSingul)
  IF(bSingul) CALL Stop_("SINGUL IN Spl_Stoikio_Calc")
  !-> inverse of A gives Elements 2:nEl in terms of components 1:nCp
  !
  ! obtain Phases 1:nFs expressed in terms of Components 1:nCp
  ALLOCATE(Y(1:nCp))
  DO I=1,nFs
    Y= tFormula(I,2:nEl)
    CALL LU_BakSub(A,Indx,Y)
    tStoikio(I,1:nCp)=Y
  ENDDO
  !
  DEALLOCATE(A)
  DEALLOCATE(Indx)
  DEALLOCATE(Y)
  !
  IF(iDebug>2) CALL Spl_Stoikio_Show
  !
  DEALLOCATE(tStoik)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Spl_Stoikio_Calc"
  !
CONTAINS

SUBROUTINE Spl_Stoikio_Show
!--
!-- write stoikio on stokio "trace file"
!--
  USE M_IoTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut,NamFInn
  !
  INTEGER:: I,iEl
  INTEGER:: FF
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< ShowStoikio"
  !
  IF(iDebug>0) WRITE(fTrc,'(A)') &
  & "Spl_Stoikio_Show -> "//TRIM(DirOut)//"_spl_stoikio.log"
  !
  CALL GetUnit(FF)
  !
  OPEN(FF,FILE=TRIM(DirOut)//"_spl_stoikio.log")
  WRITE(FF,'(A)') "!! system stoichiometry"
  !
  WRITE(FF,'(/,A,/)') "Components in terms of elements"
  DO iEl=1,nEl
    WRITE(FF,'(A1,A2,A1)',ADVANCE='NO') " ",vEle(iEl)%NamEl,T_
  ENDDO
  WRITE(FF,*)
  DO I=1,nCp
    DO iEl=1,nEl; WRITE(FF,'(F7.2,A1)',ADVANCE='NO') tStoik(I,iEl),T_; ENDDO
    WRITE(FF,'(A,A1,F15.3)') TRIM(vCpn(I)%NamCp),T_,vCpn(I)%Mole
  ENDDO
  !
  WRITE(FF,'(/,A,/)') "Phases in terms of elements"
  DO I=1,nFs
    DO iEl=1,nEl
       WRITE(FF,'(F7.3,A1)',ADVANCE='NO') tFormula(I,iEl),T_
    ENDDO
    WRITE(FF,'(A,A1,F15.3)') TRIM(vFas(I)%NamFs),T_,vFas(I)%Grt
  ENDDO
  !
  CLOSE(FF)
  !
  !CALL GetUnit(FF)
  OPEN(FF,FILE=TRIM(DirOut)//"_spl_stoikio2.log")
  WRITE(FF,'(2A)') &
  & "!! Stoikio Table of Fases (line) in terms of Components (column)", &
  & " and vFas(I)%Grt"
  !
  DO iEl=1,nCp
    WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vCpn(iEl)%NamCp),T_
  ENDDO
  WRITE(FF,*)
  
  DO I=1,nFs
    DO iEl=1,nCp
      WRITE(FF,'(F7.3,A1)',ADVANCE='NO') tStoikio(I,iEl),T_
    ENDDO
    WRITE(FF,'(A,A1,F15.3)') TRIM(vFas(I)%NamFs),T_,vFas(I)%Grt
  ENDDO
  !
  CLOSE(FF)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ ShowStoikio"
ENDSUBROUTINE Spl_Stoikio_Show

ENDSUBROUTINE Spl_Stoikio_Calc

SUBROUTINE Table_Transform(tInput,tTransform,tOutput)
!--
!-- transform, using tTransform, tInput to tOutput
!-- tTransform is invertible, (nC,nC)
!-- tInput and tOutput are matrices of same shape (nC,nS), nS>nC
!-- tOutput - inv(tTransform) * tInput
!--
  USE M_Trace,ONLY: Stop_
  USE M_Numeric_Mat,ONLY: LU_Decomp, LU_BakSub
  REAL(dp),INTENT(INOUT) :: tTransform(:,:) !square, non singular, (nC,nC)
  REAL(dp),INTENT(IN)    :: tInput(:,:)  !(nC,nS)
  REAL(dp),INTENT(OUT)   :: tOutput(:,:) !(nC,nS)
  !
  REAL(dp),ALLOCATABLE:: vY(:)
  INTEGER, ALLOCATABLE:: vIndx(:)
  LOGICAL :: bSingul
  REAL(dp):: D
  INTEGER :: I,nC,nS
  !
  nC= SIZE(tInput,1)
  nS= SIZE(tInput,2)
  !
  ALLOCATE(vY(1:nC))
  ALLOCATE(vIndx(1:nC))
  !
  !the transformation matrix, tTransform, must be invertible !!
  CALL LU_Decomp(tTransform,vIndx,D,bSingul)
  IF(bSingul) CALL Stop_("SINGUL IN Table_Transform")
  DO I= 1,nS
    vY= tInput(:,I)
    CALL LU_BakSub(tTransform,vIndx,vY)
    tOutput(:,I)= vY(:)
  ENDDO
  !
  DEALLOCATE(vY,vIndx)
  !
ENDSUBROUTINE Table_Transform

SUBROUTINE WriteRMat(f,Mat) !,N1,N2,M1,M2)
  INTEGER, INTENT(IN)::f !=File Index
  REAL(dp),INTENT(IN)::Mat(:,:)
  !
  INTEGER::I,J
  !
  DO I=1,SIZE(Mat,1)
    DO J=1,SIZE(Mat,2)
      WRITE(f,'( F7.1, A1)',ADVANCE='NO') Mat(I,J), T_
    ENDDO
    WRITE(f,*)
  ENDDO
  WRITE(f,*)
  !
ENDSUBROUTINE WriteRMat

SUBROUTINE WriteMat_I(f,Mat)
  INTEGER, INTENT(IN)::f !=File Index
  INTEGER,INTENT(IN)::Mat(:,:)
  !
  INTEGER::I,J
  !
  DO I=1,SIZE(Mat,1)
    DO J=1,SIZE(Mat,2)
      WRITE(f,'(I3, A1)',ADVANCE='NO') Mat(I,J), T_
    ENDDO
    WRITE(f,*)
  ENDDO
  WRITE(f,*)
  !
ENDSUBROUTINE WriteMat_I

SUBROUTINE Spl_WriteSystem(vCpn,tStoikio)
!--
!-- write system stoikiometry
!--
  USE M_IoTools,    ONLY: GetUnit
  USE M_Files,      ONLY: DirOut,NamFInn
  USE M_T_Component,ONLY: T_Component
  !
  USE M_Global_Vars,ONLY: tFormula,vFas
  !
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  REAL(dp),         INTENT(IN):: tStoikio(:,:)

  INTEGER :: nCp,nFs,iCpn,iFs
  INTEGER :: FF
  REAL(dp):: X
  !
  nCp=SIZE(vCpn)
  nFs=SIZE(vFas)
  !
  WRITE(fTrc,'(/,A)') "< WriteSystem"
  WRITE(fTrc,'(A)') "Spl_WriteSystem -> "//TRIM(DirOut)//"_spl_system.log"

  CALL GetUnit(FF)
  OPEN(FF,FILE=TRIM(DirOut)//"_spl_system.log")

  WRITE(FF,'(A,A1)',ADVANCE='NO') "_", T_
  DO iCpn=1,nCp
    WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vCpn(iCpn)%NamCp), T_
  ENDDO
  WRITE(FF,'(A)') "GIBBS"

  WRITE(FF,'(A,A1)',ADVANCE='NO') "MOLES", T_
  DO iCpn=1,nCp
    WRITE(FF,'(G12.3,A1)',ADVANCE='NO') vCpn(iCpn)%Mole, T_
  ENDDO
  WRITE(FF,'(A)') "0.000"

  DO iFs=1,nFs
    WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vFas(iFs)%NamFs), T_
    DO iCpn=1,nCp
      WRITE(FF,'(F7.2,A1)',ADVANCE='NO') tStoikio(iFs,iCpn), T_
    ENDDO
    X= vFas(iFs)%Grt
    IF(tFormula(iFs,1)/=0) X= vFas(iFs)%Grt /tFormula(iFs,1)
    ! WRITE(FF,'(2F15.3)') vFas(iFs)%Grt,X
    WRITE(FF,'(2F15.3)') vFas(iFs)%Grt
  ENDDO

  WRITE(fTrc,'(/,A,/)') "</ WriteSystem"

  CLOSE(FF)

  RETURN
END SUBROUTINE Spl_WriteSystem

ENDMODULE M_Simplex_Build

