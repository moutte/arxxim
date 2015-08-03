MODULE M_Solmodel_Pitzer_Dtb
!--
!-- Module routines associées à Pitzer
!-- Livraison 10/2006: Nicolas Ferrando
!-- Remoduling : Anthony Michel
!--
!-- modified : J.Moutte, 11/2006
!-- -> M_Pitzer_Dtb-  database reading, and TP computations
!-- -> M_Pitzer_Calc- calculations
!--
  USE M_Kinds
  USE M_Trace,    ONLY: fTrc,iDebug,T_,Stop_
  USE M_T_Species,ONLY: T_Species
  IMPLICIT NONE
  !
  PRIVATE
  !
  !-- public --
  PUBLIC:: Solmodel_Pitzer_Dtb_Init
  PUBLIC:: Solmodel_Pitzer_Dtb_TPUpdate
  PUBLIC:: Solmodel_Pitzer_Dtb_Clean
  PUBLIC:: Solmodel_Pitzer_Dtb_TPtest
  !
  REAL(dp),DIMENSION(:,:),  ALLOCATABLE,PUBLIC:: tBeta0,tBeta1,tBeta2
  REAL(dp),DIMENSION(:,:),  ALLOCATABLE,PUBLIC:: tCPhi,tTheta,tLamda
  REAL(dp),DIMENSION(:,:),  ALLOCATABLE,PUBLIC:: tAlfa1,tAlfa2
  REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,PUBLIC:: tZeta,tPsi
  !
  !-- private --
  TYPE(T_Species),ALLOCATABLE:: vSolut(:)
  !
  INTEGER,DIMENSION(:,:),  ALLOCATABLE:: tI_Beta0,tI_Beta1,tI_Beta2
  INTEGER,DIMENSION(:,:),  ALLOCATABLE:: tI_CPhi,tI_Theta,tI_Lamda
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: tI_Zeta,tI_Psi
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vBeta0,vBeta1,vBeta2
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vCPhi,vTheta,vLamda
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vZeta,vPsi
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vAlfa1,vAlfa2
  !
  ! usage:
  !  k= tI_Beta0(i,j)
  !  IF(k/=0) THEN
  !    tBeta0(i,j)= vBeta0(k)

  ! FUNCTION CoeffAtTP()
  !    Res= T%X0
  !    &  + T%X3 *(TdgK -Tref) &
  !    &  + T%X1 *(One/TdgK -One/Tref) &
  !    &  + T%x2 *log(TdgK/Tref)       &
  !    &  + T%x4 *(TdgK*TdgK -Tref*Tref)

  !~ !! former version: static fitting parameter
  !~ TYPE:: T_Fitt
    !~ INTEGER:: iModel
    !~ REAL(dp):: vX(1:8)
  !~ ENDTYPE T_Fitt

  TYPE:: T_Fitt
    CHARACTER(LEN=30):: namModel
    INTEGER:: iModel ! model selector
    INTEGER:: N      ! dimension of vX (model- dependent)
    REAL(dp),ALLOCATABLE:: vX(:)
  ENDTYPE T_Fitt
  !
  TYPE(T_Fitt),DIMENSION(:),ALLOCATABLE:: vF_Beta0,vF_Beta1,vF_Beta2
  TYPE(T_Fitt),DIMENSION(:),ALLOCATABLE:: vF_CPhi,vF_Theta,vF_Lamda
  TYPE(T_Fitt),DIMENSION(:),ALLOCATABLE:: vF_Zeta,vF_Psi

CONTAINS

SUBROUTINE Solmodel_Pitzer_Dtb_Init(vSpc)
  USE M_T_Species, ONLY: T_Species
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  !
  LOGICAL:: Error
  CHARACTER(LEN=80):: ErrorMsg
  !
  CALL Solmodel_Pitzer_Dtb_Clean
  !
  CALL Solmodel_Pitzer_Dtb_Read(vSpc,Error,ErrorMsg)
  !
  IF(Error) &
  & CALL Stop_("Error in Solmodel_PitzerDtb_Read:"//TRIM(ErrorMsg))
  !
  !CALL Solmodel_PitzerDtb_TPUpdate(TdgK)
  !
ENDSUBROUTINE Solmodel_Pitzer_Dtb_Init

SUBROUTINE Solmodel_Pitzer_Dtb_Read( &
& vSpc,           &
& Error,ErrorMsg  &
& )
!--
!-- read interaction coeffs from the Pitzer database
!-- must be called AFTER vSpc has been built
!-- -> vPitz will be a subset of vSpc
!--
  USE M_IoTools
  USE M_Files,        ONLY: NamFPtz
  USE M_Numeric_Const,ONLY: ln10
  USE M_T_Species,    ONLY: T_Species,Species_Index
  !
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  !
  LOGICAL,         INTENT(OUT):: Error
  CHARACTER(LEN=*),INTENT(OUT):: ErrorMsg
  !
  REAL(dp),PARAMETER:: Tref= 298.15D0
  !
  LOGICAL:: L_ScanSpecies
  !
  CHARACTER(LEN=512):: L,W,W2
  LOGICAL           :: EoL
  !
  ! CHARACTER(LEN=7):: FittFormat
  CHARACTER(LEN=30):: namModel
  !
  LOGICAL :: Ok
  INTEGER :: vISpc(3)
  INTEGER :: iSp1,iSp2,iSp3,nSl,I,J,N
  INTEGER :: F,ios,ff
  INTEGER :: iModel
  INTEGER :: iBeta0,iBeta1,iBeta2
  INTEGER :: iCPhi,iTheta,iLamda
  INTEGER :: iPsi,iZeta
  REAL(dp):: vXread(dimV) !,X0,vX(dimV)
  TYPE(T_Fitt):: Coeff
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Solmodel_PitzerDtb_Read"
  !
  CALL Solmodel_Pitzer_Dtb_Clean
  ! Pitzer_Dtb_Init_Done=.FALSE.
  !
  L_ScanSpecies= .FALSE.
  !
  Error= .FALSE.
  ErrorMsg= "OK"
  !
  !!  IF(PRESENT(TdgK)) THEN; T=TdgK
  !!  ELSE;                   T=298.15D0
  !!  ENDIF
  !!  IF(PRESENT(Pbar)) THEN; P=Pbar
  !!  ELSE;                   P=1.01325D0 !1 atm in bars
  !!                          !caveat: some FUNCTIONs may USE Pascal ?...
  !!  ENDIF
  !!
  !! nSl= SIZE(vSpc)
  !
  ! nSl is number of solute species
  nSl= COUNT(vSpc(:)%Typ=="AQU") - 1 !mod 19/06/2008 09:38
  !
  ALLOCATE(vSolut(1:nSl))
  !
  J= 0
  DO I=1,SIZE(vSpc)
    IF(vSpc(I)%Typ=="AQU" .AND. vSpc(I)%NamSp/="H2O") THEN
      J= J+1
      vSolut(J)= vSpc(I)
    ENDIF
  ENDDO
  !
  iModel= 0
  namModel= "NONE"
  !
  ALLOCATE(tBeta0(1:nSl,1:nSl))      ;  tBeta0=Zero
  ALLOCATE(tBeta1(1:nSl,1:nSl))      ;  tBeta1=Zero
  ALLOCATE(tBeta2(1:nSl,1:nSl))      ;  tBeta2=Zero
  ALLOCATE(tCPhi (1:nSl,1:nSl))      ;  tCPhi= Zero
  ALLOCATE(tTheta(1:nSl,1:nSl))      ;  tTheta=Zero
  ALLOCATE(tLamda(1:nSl,1:nSl))      ;  tLamda=Zero
  ALLOCATE(tAlfa1(1:nSl,1:nSl))      ;  tAlfa1=Zero
  ALLOCATE(tAlfa2(1:nSl,1:nSl))      ;  tAlfa2=Zero
  ALLOCATE(tZeta (1:nSl,1:nSl,1:nSl));  tZeta= Zero
  ALLOCATE(tPsi  (1:nSl,1:nSl,1:nSl));  tPsi=  Zero
  !
  ALLOCATE(tI_Beta0(1:nSl,1:nSl))    ;  tI_Beta0= 0
  ALLOCATE(tI_Beta1(1:nSl,1:nSl))    ;  tI_Beta1= 0
  ALLOCATE(tI_Beta2(1:nSl,1:nSl))    ;  tI_Beta2= 0
  ALLOCATE(tI_CPhi (1:nSl,1:nSl))    ;  tI_CPhi=  0
  ALLOCATE(tI_Theta(1:nSl,1:nSl))    ;  tI_Theta= 0
  ALLOCATE(tI_Lamda(1:nSl,1:nSl))    ;  tI_Lamda= 0

  !ALLOCATE(tI_Alfa1(1:nSl,1:nSl));        tI_Alfa1=.false.
  !ALLOCATE(tI_Alfa2(1:nSl,1:nSl));        tI_Alfa2=.false.

  ALLOCATE(tI_Zeta (1:nSl,1:nSl,1:nSl)) ; tI_Zeta= 0
  ALLOCATE(tI_Psi  (1:nSl,1:nSl,1:nSl)) ; tI_Psi=  0
  !
  CALL GetUnit(F)
  !
  !----- reading database, FIRST PASS -> read tI_Beta0, tI_Beta1, ... --
  OPEN(F,FILE=TRIM(NamFPtz),STATUS='OLD')

  DoFile1: DO

    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile1
    !
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE doFile1 !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT DoFile1
    !
    CASE("PITZER")
      !
      doLine1: DO
        !
        READ(F,'(A)',IOSTAT=ios) L
        IF(ios/=0) EXIT doFile1
        !
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE doLine1 !skip comment lines
        CALL AppendToEnd(L,W,EoL)

        SELECT CASE(W)

        CASE("ENDINPUT");        EXIT doFile1
        CASE("ENDPITZER","END"); EXIT doLine1

        CASE("FORMAT")
          CALL LinToWrd(L,W2,EoL)
          SELECT CASE(TRIM(W2))
          CASE("NEW")  ;  L_ScanSpecies= .TRUE.
          CASE("OLD")  ;  L_ScanSpecies= .FALSE.
          CASE DEFAULT
            Error= .TRUE.
            ErrorMsg= TRIM(W2)//"= not valid as FORMAT"
            RETURN
          END SELECT

        CASE("FITTING")
          !CALL LinToWrd(L,W2,EoL)
          CYCLE doLine1
        !
        !--------------------------------------- 2-species parameters --
        CASE("BETA0","BETA1", "BETA2","CPHI", "LAMBDA","THETA")
          IF(L_ScanSpecies) THEN
            CALL LinToWrd(L,W2,EoL)
            CALL Species_Scan(W2,2,vISpc,Error)
            IF(Error) THEN
              ErrorMsg= TRIM(W2)//" -> Error in species list ?"
              RETURN
            ELSE
              iSp1= vISpc(1)
              iSp2= vISpc(2)
            ENDIF
          ELSE
            CALL LinToWrd(L,W2,EoL)  ;  iSp1= Species_Index(W2,vSolut)
            CALL LinToWrd(L,W2,EoL)  ;  iSp2= Species_Index(W2,vSolut)
          ENDIF
          !
          !--- skip the lines with species not in current species set --
          IF(iSp1*iSp2==0) CYCLE doLine1
          !---/
          !
          IF(TRIM(W)=="LAMBDA") THEN
            IF(vSolut(iSp2)%Z/=0) THEN
              ErrorMsg= "LAMBDA coeff: species 2 must be neutral"
              IF(iDebug>2) WRITE(fTrc,'(A,/,A)'), &
              & TRIM(ErrorMsg), &
              & TRIM(vSolut(iSp1)%NamSp)//"="// &
              & TRIM(vSolut(iSp2)%NamSp)
              Error= .TRUE.
              RETURN
            END IF
          ELSE
            CALL Sort2sp(iSp1,iSp2,Error)
          END IF
          !
          IF(Error) &
          & CALL Stop_("Pitzer database: error in line "//TRIM(W)//TRIM(W2))
          !
          !IF(iSp1*iSp2>1) THEN !-> the 2 species are in database
          SELECT CASE(TRIM(W))
            CASE("BETA0")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Beta0)
            CASE("BETA1")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Beta1)
            CASE("BETA2")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Beta2)
            CASE("CPHI")      ; Ok= Ok_2Sp(iSp1,iSp2,tI_CPhi )
            CASE("THETA")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Theta)
            CASE("LAMBDA")    ; Ok= Ok_2Sp(iSp1,iSp2,tI_Lamda)
            !! CASE("Alfa") !! not used
            !!! Alfa is calculated according to Harvie-etal-84
            !!! see also Pitzer, ReviewInMineralogy-17, p104
          END SELECT
            !
          !ENDIF
        !ENDCASE Beta...
        !--------------------------------------/ 2-species parameters --
        !
        !--------------------------------------- 3-species parameters --
        CASE("PSI","ZETA")
          !
          IF(L_ScanSpecies) THEN
            CALL LinToWrd(L,W2,EoL)
            CALL Species_Scan(W2,3,vISpc,Error)
            IF(Error) THEN
              ErrorMsg= TRIM(W2)//" -> Error in species list ?"
              RETURN
            ELSE
              iSp1= vISpc(1)
              iSp2= vISpc(2)
              iSp3= vISpc(3)
            ENDIF
          ELSE
            CALL LinToWrd(L,W2,EoL); iSp1=Species_Index(W2,vSolut)
            CALL LinToWrd(L,W2,EoL); iSp2=Species_Index(W2,vSolut)
            CALL LinToWrd(L,W2,EoL); iSp3=Species_Index(W2,vSolut)
          ENDIF
          !
          IF(iSp1*iSp2*iSp3==0) CYCLE doLine1
          !
          IF(TRIM(W)=="PSI") THEN
            IF(vSolut(iSp1)%Z*vSolut(iSp2)%Z<0) THEN
              ErrorMsg= "PSI coeff: Species 1 & 2 must be ions same charge"
              IF(iDebug>2) WRITE(fTrc,'(A,/,A)'), &
              & TRIM(ErrorMsg), &
              & TRIM(vSolut(iSp1)%NamSp)//"="// &
              & TRIM(vSolut(iSp2)%NamSp)//"="// &
              & TRIM(vSolut(iSp3)%NamSp)
              Error= .TRUE.
              RETURN
            END IF
          END IF
          !
          IF(TRIM(W)=="ZETA") THEN
            IF(vSolut(iSp1)%Z*vSolut(iSp2)%Z>0) THEN
              ErrorMsg= "PSI coeff: Species 1 & 2 must be ions opposite charge"
              IF(iDebug>2) WRITE(fTrc,'(A,/,A)'), &
              & TRIM(ErrorMsg), &
              & TRIM(vSolut(iSp1)%NamSp)//"="// &
              & TRIM(vSolut(iSp2)%NamSp)//"="// &
              & TRIM(vSolut(iSp3)%NamSp)
              Error= .TRUE.
              RETURN
            END IF
          END IF
          !
          ! CALL Sort3sp(iSp1,iSp2,iSp3,Error)
          !
          CALL Sort2sp(iSp1,iSp2,Error)
          !
          IF(Error) &
          & CALL Stop_("Pitzer database: error in line "//TRIM(W)//TRIM(W2))
          !
          SELECT CASE(W)
            !CASE("PSI");   tI_Psi (iSp1,iSp2,iSp3)= 1
            !CASE("ZETA");  tI_Zeta(iSp1,iSp2,iSp3)= 1
            CASE("PSI")  ;  Ok= Ok_3Sp(iSp1,iSp2,iSp3,tI_Psi )
            CASE("ZETA") ;  Ok= Ok_3Sp(iSp1,iSp2,iSp3,tI_Zeta)
          END SELECT
          !
        !ENDCASE PSI...
        !--------------------------------------/ 3-species parameters --
        CASE DEFAULT
          Error= .TRUE.
          ErrorMsg= TRIM(W)//" is unknown keyword"
          RETURN
        !
        END SELECT

      ENDDO doLine1
    !END CASE("PITZER")
    
    END SELECT

  ENDDO DoFile1
  !
  CLOSE(F)
  !----/ reading database, FIRST PASS -> read tI_Beta0, tI_Beta1, ... --

  !------------------------------------------------------ allocations --
  iBeta0= COUNT(tI_Beta0(:,:)==1)
  iBeta1= COUNT(tI_Beta1(:,:)==1)
  iBeta2= COUNT(tI_Beta2(:,:)==1)
  iCPhi=  COUNT(tI_CPhi (:,:)==1)
  iTheta= COUNT(tI_Theta(:,:)==1)
  iLamda= COUNT(tI_Lamda(:,:)==1)
  !
  ALLOCATE(vBeta0(iBeta0))  ;  ALLOCATE(vF_Beta0(iBeta0))
  ALLOCATE(vBeta1(iBeta1))  ;  ALLOCATE(vF_Beta1(iBeta1))
  ALLOCATE(vBeta2(iBeta2))  ;  ALLOCATE(vF_Beta2(iBeta2))
  !
  ALLOCATE(vAlfa1(iBeta1))
  ALLOCATE(vAlfa2(iBeta2))
  !
  ALLOCATE(vCPhi (iCPhi) )  ;  ALLOCATE(vF_CPhi (iCPhi) )
  ALLOCATE(vTheta(iTheta))  ;  ALLOCATE(vF_Theta(iTheta))
  ALLOCATE(vLamda(iLamda))  ;  ALLOCATE(vF_Lamda(iLamda))
  !
  iBeta0= 0
  iBeta1= 0
  iBeta2= 0
  iCPhi=  0
  iTheta= 0
  iLamda= 0
  !
  iPsi=  COUNT(tI_Psi(:,:,:)==1)
  iZeta= COUNT(tI_Zeta(:,:,:)==1)
  ALLOCATE(vPsi(iPsi))     ;  ALLOCATE(vF_Psi(iPsi))
  ALLOCATE(vZeta(iZeta))   ;  ALLOCATE(vF_Zeta(iZeta))
  !
  iPsi=  0
  iZeta= 0
  !------------------------------------------------------/allocations --

  !-------- reading database, SECOND PASS -> read vBeta0, vBeta1, ... --
  OPEN(F,FILE=TRIM(NamFPtz),STATUS='OLD')
  !
  DoFile2: DO

    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT doFile2
    !
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE doFile2 !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT doFile2
    !
    CASE("PITZER")
      !
      doLine2: DO
        !
        READ(F,'(A)',IOSTAT=ios) L
        IF(ios/=0) EXIT doFile2
        !
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE doLine2 !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT");        EXIT doFile2
        CASE("ENDPITZER","END"); EXIT doLine2
        !
        CASE("FITTING")
          CALL LinToWrd(L,W2,EoL)
          namModel= TRIM(W2)
          SELECT CASE(TRIM(W2))
            CASE("NONE")               ;  iModel= 0
            CASE("PHREEQC")            ;  iModel= 1
            CASE("KUEHN")              ;  iModel= 2
            CASE("CHRISTOV-2004")      ;  iModel= 3
            CASE("PITZER-1984")        ;  iModel= 4
            CASE("PABALAN-1987")       ;  iModel= 5
            CASE("POLYA-2007")         ;  iModel= 6
            CASE("LI-DUAN-2007")       ;  iModel= 7
            CASE DEFAULT
              Error= .TRUE.
              ErrorMsg= TRIM(W2)//"= unknown keyword for FITTING"
              RETURN
          END SELECT
        !
        !------------------------------------- 2-species parameters --
        CASE("BETA0","BETA1", "BETA2","CPHI", "LAMBDA","THETA")
          !------------------------------- retrieve species indexes --
          IF(L_ScanSpecies) THEN
            CALL LinToWrd(L,W2,EoL)
            CALL Species_Scan(W2,2,vISpc,Error)
            IF(Error) THEN
              ErrorMsg= TRIM(W2)//" -> Error in species list ?"
            ELSE
              iSp1= vISpc(1)
              iSp2= vISpc(2)
            ENDIF
          ELSE
            CALL LinToWrd(L,W2,EoL)  ;  iSp1= Species_Index(W2,vSolut)
            CALL LinToWrd(L,W2,EoL)  ;  iSp2= Species_Index(W2,vSolut)
          ENDIF
          !
          IF(iSp1*iSp2==0) CYCLE doLine2
          !
          CALL Sort2sp(iSp1,iSp2,Error)
          !--------------------------------/ retrieve species indexes --
          !
        !--------------------------------------- 3-species parameters --
        CASE("PSI","ZETA")
          !--------------------------------- retrieve species indexes --
          IF(L_ScanSpecies) THEN
            CALL LinToWrd(L,W2,EoL)
            CALL Species_Scan(W2,3,vISpc,Error)
            IF(Error) THEN
              ErrorMsg= TRIM(W2)//" -> Error in species list ?"
              RETURN
            ELSE
              iSp1= vISpc(1)
              iSp2= vISpc(2)
              iSp3= vISpc(3)
            ENDIF
          ELSE
            CALL LinToWrd(L,W2,EoL); iSp1=Species_Index(W2,vSolut)
            CALL LinToWrd(L,W2,EoL); iSp2=Species_Index(W2,vSolut)
            CALL LinToWrd(L,W2,EoL); iSp3=Species_Index(W2,vSolut)
          ENDIF
          !
          IF(iSp1*iSp2*iSp3==0) CYCLE doLine2
          !
          !! CALL Sort3sp(iSp1,iSp2,iSp3,Error)
          CALL Sort2sp(iSp1,iSp2,Error)
          !
          !--------------------------------/ retrieve species indexes --
        CASE DEFAULT
          Error= .TRUE.
          ErrorMsg= TRIM(W)//" is unknown keyword"
        !
        END SELECT
        !
        !-------------------------------------------- read parameters --
        !
        CALL ReadRValsV(L,N,vXread)
        !
        SELECT CASE(iModel)
          CASE(0)  ;  N= 1  ! "NONE"
          CASE(1)  ;  N= 5  ! "PHREEQ"
          CASE(2)  ;  N= 5  ! "KUEHN"
          CASE(3)  ;  N= 9  ! "CHRISTOV-2004"
          CASE(4)  ;  N= 21 ! "PITZER-1984", Na+=Cl-
          CASE(5)  ;  N= 12 ! "PABALAN-1987"
          CASE(6)  ;  N= 8  ! "POLYA-2007"
          CASE(7)  ;  N= 11 ! "LI-DUAN-2007"
          CASE DEFAULT
            CALL Stop_("Error on Pitz_Dtb / iModel")
        END SELECT
        !
        ALLOCATE(Coeff%vX(N))
        !
        !!--- for static fitting model
        !Coeff%vX(1:8)= vXread(1:8)
        !Coeff%iModel= iModel
        !!---/
        Coeff%N= N
        Coeff%namModel= TRIM(namModel)
        Coeff%iModel= iModel
        Coeff%vX(1:N)= vXread(1:N)
        !
        SELECT CASE(W)
      
        !--------------------------------------- 2-species parameters --
        CASE("BETA0")
          iBeta0= iBeta0 +1
          tI_Beta0(iSp1,iSp2)= iBeta0
          ALLOCATE(vF_Beta0(iBeta0)%vX(Coeff%N))
          vF_Beta0(iBeta0)= Coeff
        CASE("BETA1")
          iBeta1= iBeta1 +1
          tI_Beta1(iSp1,iSp2)= iBeta1
          ALLOCATE(vF_Beta1(iBeta1)%vX(Coeff%N))
          vF_Beta1(iBeta1)= Coeff
        CASE("BETA2")
          iBeta2= iBeta2 +1
          tI_Beta2(iSp1,iSp2)= iBeta2
          ALLOCATE(vF_Beta2(iBeta2)%vX(Coeff%N))
          vF_Beta2(iBeta2)= Coeff
        CASE("CPHI")
          iCPhi= iCPhi +1
          tI_CPhi(iSp1,iSp2)= iCPhi
          ALLOCATE(vF_CPhi(iCPhi)%vX(Coeff%N))
          vF_CPhi(iCPhi)= Coeff
        CASE("THETA")
          iTheta= iTheta +1
          tI_Theta(iSp1,iSp2)= iTheta
          ALLOCATE(vF_Theta(iTheta)%vX(Coeff%N))
          vF_Theta(iTheta)= Coeff
        CASE("LAMBDA")
          iLamda= iLamda +1
          tI_Lamda(iSp1,iSp2)= iLamda
          ALLOCATE(vF_Lamda(iLamda)%vX(Coeff%N))
          vF_Lamda(iLamda)= Coeff

        !! CASE("Alfa") !! not used
        !!! Alfa is calculated according to Harvie-etal-84
        !!! see also Pitzer, ReviewInMineralogy-17, p104

        !--------------------------------------- 3-species parameters --
        CASE("PSI")
          iPsi= iPsi+1
          !vPsi(iPsi)= X
          tI_Psi(iSp1,iSp2,iSp3)= iPsi
          ALLOCATE(vF_Psi(iPsi)%vX(Coeff%N))
          vF_Psi(iPsi)= Coeff
        CASE("ZETA")
          iZeta= iZeta +1
          tI_Zeta(iSp1,iSp2,iSp3)= iZeta
          ALLOCATE(vF_Zeta(iZeta)%vX(Coeff%N))
          vF_Zeta(iZeta)= Coeff

        END SELECT
        !
        DEALLOCATE(Coeff%vX)
        !-------------------------------------------/ read parameters --
        !
      ENDDO doLine2
    !ENDCASE("BLOCK")
    
    END SELECT

  ENDDO doFile2
  !
  CLOSE(F)
  !--------/reading database, SECOND PASS -> read vBeta0, vBeta1, ... --
  !
  !----------------------------------------------------- Alfa1,tAlfa2 --
  !--- following Pitzer-RevMin17-p104, Harvie-etal-84
  !--- cf CHRISTOV-2004-Moeller-GCA68-3718
  DO iSp1=1,nSl
    DO iSp2=1,nSl
      IF(tI_Beta1(iSp1,iSp2)>0) THEN
        IF( ABS(vSolut(iSp1)%Z)<2 .OR. ABS(vSolut(iSp2)%Z)<2 ) THEN
          vAlfa1(tI_Beta1(iSp1,iSp2))= 2.0D0
        ELSE
          vAlfa1(tI_Beta1(iSp1,iSp2))= 1.4D0
        ENDIF
      ENDIF

      IF(tI_Beta2(iSp1,iSp2)>0) THEN
        IF(ABS(vSolut(iSp1)%Z)<2 .OR. ABS(vSolut(iSp2)%Z)<2) THEN
          !-> ignore the tBeta2 term in B_MX expression
          vAlfa2(tI_Beta2(iSp1,iSp2))= Zero
        ELSE
          vAlfa2(tI_Beta2(iSp1,iSp2))= 12.0D0
        ENDIF
      ENDIF
      !
    ENDDO
  ENDDO
  !----------------------------------------------------/ Alfa1,tAlfa2 --
  !
  CALL Solmodel_Pitzer_Dtb_TPUpdate(298.15D0,1.0D0)
  !
  IF(iDebug>0) THEN !=========================================< trace ==
    CALL GetUnit(ff)
    OPEN(ff,FILE="debug_pitzer_data.log")
    !
    CALL Solmodel_Pitzer_Dtb_Check(ff,vSolut)
    CALL Solmodel_Pitzer_Dtb_Check_Old(ff,vSolut)
    !
    CLOSE(ff)
  ENDIF !=====================================================< trace ==
  ! Solmodel_Pitzer
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Solmodel_Pitzer_ Solmodel_PitzerDtb_Read"
  !
CONTAINS

LOGICAL FUNCTION Ok_2Sp(I1,I2,V)
!--
!-- check that a species pair is not already read,
!-- update the table V
!--
  INTEGER,INTENT(IN)   :: I1,I2
  INTEGER,INTENT(INOUT):: V(:,:)
  !
  Ok_2Sp= (V(I1,I2)==0)
  IF(Ok_2Sp) THEN
    V(I1,I2)= 1
  ELSE
    PRINT '(2(A,1X))',TRIM(vSolut(I1)%NamSp),TRIM(vSolut(I2)%NamSp)
    CALL Stop_("Pitz_Dtb.Error: already allocated")
  ENDIF
  !
END FUNCTION Ok_2Sp

LOGICAL FUNCTION Ok_3Sp(I1,I2,I3,V)
!--
!-- check that a species triplet is not already read,
!-- update the table V
!--
  INTEGER,INTENT(IN)   :: I1,I2,I3
  INTEGER,INTENT(INOUT):: V(:,:,:)
  !
  Ok_3Sp= (V(I1,I2,I3)==0)
  IF(Ok_3Sp) THEN
    V(I1,I2,I3)= 1
  ELSE
    CALL Stop_("Error in Pitz_Dtb: triplet already allocated")
  ENDIF
  !
END FUNCTION Ok_3Sp

END SUBROUTINE Solmodel_Pitzer_Dtb_Read

SUBROUTINE Sort2sp(I1,I2,Error)
!--
!-- order the I1,I2 pair,
!-- to put the data in the upper triangle
!--
  INTEGER,INTENT(INOUT):: I1,I2
  LOGICAL,INTENT(OUT)  :: Error
  !
  INTEGER:: J1,J2
  !
  Error= (I1==I2)
  !
  IF(Error) RETURN
  !
  J1= MIN(I1,I2)  ;  J2= MAX(I1,I2)
  I1= J1          ;  I2= J2
  !
END SUBROUTINE Sort2sp

SUBROUTINE Sort3sp(I1,I2,I3,Error)
!-- order I1,I2,I3 to I1-I2-I3
  INTEGER,INTENT(INOUT):: I1,I2,I3
  LOGICAL,INTENT(OUT)  :: Error
  !
  INTEGER:: J1,J2,J3
  !
  Error= (I1==I2 .OR. I1==I3 .OR. I2==I3)
  IF(Error) RETURN
  !
  J1= MIN(I1,I2,I3)
  J3= MAX(I1,I2,I3)
  !
  IF(I1/=J1 .AND. I1/=J3) J2= I1
  IF(I2/=J1 .AND. I2/=J3) J2= I2
  IF(I3/=J1 .AND. I3/=J3) J2= I3
  !
  I1= J1  ;  I2= J2  ;  I3= J3
  !
END SUBROUTINE Sort3sp

SUBROUTINE Species_Scan(Str,N,iSpc,Error)
!--
!-- scan a species string, where species are separated by "-",
!-- and compute corresponding indexes in vSolut
!--
  USE M_T_Species,ONLY: Species_Index
  CHARACTER(LEN=*),INTENT(IN):: Str
  INTEGER,INTENT(IN) :: N
  INTEGER,INTENT(OUT):: iSpc(:)
  LOGICAL,INTENT(OUT):: Error
  !
  CHARACTER(LEN=23):: Str1,Str2
  INTEGER:: I,J,K
  !
  Error= .FALSE.
  !
  Str1= TRIM(Str)//"=" !append a separator at end
  Str2= ""
  K= 0
  DO
    !WRITE(91,'(2A)') "str1=",TRIM(Str1)
    I= INDEX(Str1,'=') ! first occurence of separator
    J= LEN_TRIM(Str1)
    IF(I==0) EXIT
    WRITE(Str2,'(A)') TRIM(str1(1:I-1)) ! Str2= left of first separator
    WRITE(Str1,'(A)') TRIM(str1(I+1:J)) ! Str1= right of first separator
    K= K+1
    iSpc(K)= Species_Index(Str2,vSolut)
    !WRITE(91,'(A,I3,1X,A)') "    I,str2=",iSpc(K),TRIM(Str2)
    IF(K==N) EXIT
  ENDDO
  !pause
  !WRITE(91,*)
  !
  Error= (K<N)
  !
  RETURN
END SUBROUTINE Species_Scan

SUBROUTINE Solmodel_Pitzer_Dtb_TPUpdate(TdgK,Pbar)
  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Beta0, vBeta0)
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Beta1, vBeta1)
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Beta2, vBeta2)
  !
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_CPhi,  vCPhi)
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Theta, vTheta)
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Lamda, vLamda)
  !
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Psi,  vPsi)
  CALL Fitt_TPUpdate(TdgK,Pbar, vF_Zeta, vZeta)
  !
  CALL Pitzer_TableFill
  !
  RETURN
END SUBROUTINE  Solmodel_Pitzer_Dtb_TPUpdate

SUBROUTINE Pitzer_TableFill
  !
  CALL TableFill(tI_Beta0, vBeta0, tBeta0)
  CALL TableFill(tI_Beta1, vBeta1, tBeta1)
  CALL TableFill(tI_Beta2, vBeta2, tBeta2)
  !
  CALL TableFill(tI_Beta1, vAlfa1, tAlfa1)
  CALL TableFill(tI_Beta2, vAlfa2, tAlfa2)
  !
  CALL TableFill(tI_CPhi,  vCPhi,  tCPhi )
  CALL TableFill(tI_Theta, vTheta, tTheta)
  CALL TableFill(tI_Lamda, vLamda, tLamda)
  !
  CALL TableFill3(tI_Psi,  vPsi,  tPsi)
  CALL TableFill3(tI_Zeta, vZeta, tZeta)
  !
  RETURN
END SUBROUTINE Pitzer_TableFill

SUBROUTINE TableFill(tI_In,vCoef,tR_Out)
  INTEGER, INTENT(IN) :: tI_In(:,:)
  REAL(dp),INTENT(IN) :: vCoef(:)
  REAL(dp),INTENT(OUT):: tR_Out(:,:)
  !
  INTEGER:: I,J,K
  !
  DO I=1,SIZE(tI_In,1)
    DO J=1,SIZE(tI_In,2)
      K= tI_In(I,J)
      IF(K>0) tR_Out(I,J)= vCoef(K)
    ENDDO
  ENDDO
  !
END SUBROUTINE TableFill

SUBROUTINE TableFill3(tI_In,vCoef,tR_Out)
  INTEGER, INTENT(IN) :: tI_In(:,:,:)
  REAL(dp),INTENT(IN) :: vCoef(:)
  REAL(dp),INTENT(OUT):: tR_Out(:,:,:)
  !
  INTEGER:: I,J,K,L
  !
  DO I=1,SIZE(tI_In,1)
    DO J=1,SIZE(tI_In,2)
      DO K=1,SIZE(tI_In,3)
        L= tI_In(I,J,K)
        IF(L>0) tR_Out(I,J,K)= vCoef(L)
      ENDDO
    ENDDO
  ENDDO
  !
END SUBROUTINE TableFill3

SUBROUTINE Fitt_TPUpdate(T,P,vFitt,vCoef)
  REAL(dp),INTENT(IN):: T,P !TdgK,Pbar
  TYPE(T_Fitt),INTENT(IN):: vFitt(:)
  REAL(dp),INTENT(OUT):: vCoef(:)
  !
  INTEGER:: I
  REAL(dp),PARAMETER:: Tref= 298.15D0
  TYPE(T_Fitt):: Fit
  !
  DO I=1,SIZE(vCoef)
    !
    Fit= vFitt(I)
    !
    SELECT CASE(Fit%iModel)

    CASE(0) !"NONE"
      vCoef(I)= Fit%vX(1) !"NONE"

    CASE(1) !"PHREEQ"
      vCoef(I)= Fit%vX(1) &
      &       + Fit%vX(2) *(One/T -One/Tref) &
      &       + Fit%vX(3) *LOG(T/Tref)       &
      &       + Fit%vX(4) *(T     -Tref)     &
      &       + Fit%vX(5) *(T*T   -Tref*Tref)
      !
    CASE(2) !"KUEHN"
      vCoef(I)= Fit%vX(1) &
      &       + Fit%vX(2) *(One/T -One/Tref) &
      &       + Fit%vX(3) *LOG(T/Tref)       &
      &       + Fit%vX(4) *(T     -Tref)     &
      &       + Fit%vX(5) *(T*T   -Tref*Tref)
      !
    CASE(3) !"CHRISTOV-2004"
      vCoef(I)= Fit%vX(1) &              ! a1
      &       + Fit%vX(2) *T           & ! a2*T
      &       + Fit%vX(3) *T*T         & ! a3*T**2
      &       + Fit%vX(4) *T**3        & ! a4*T**3
      &       + Fit%vX(5) *(One/T)     & ! a5/T
      &       + Fit%vX(6) *LOG(T)      & ! a6*ln(T)
      &       + Fit%vX(7) /(T-263.0D0) & ! a7/(T-263)
      &       + Fit%vX(8) /(680.0D0-T) & ! a8/(680-T)
      &       + Fit%vX(9) /(T-227.0D0)   ! a9/(T-227)
      !
      ! AZAROUAL, Chemical Geology,
      ! f(T)= a1 + a2 T  + a3 /T         + a4 lnT       + a5 /(T- 263)
      !          + a6 T² + a7 /(680 - T) + a8 /(T- 227)
      ! nearly the same equation as Christov, 2004,
      ! with a different order and an additional term in /(T- 227)
      ! a1 = 1
      ! a2 = 2
      ! a3 = 5
      ! a4 = 6
      ! a5 = 7
      ! a6 = 3
      ! a7 = 8
      ! a8 = 9
      !
    CASE(4) !"PITZER-1984", Na+=Cl-
      ! K.S. Pitzer, J.C. Peiper and R.H. Busey, 1984
      ! Thermodynamic properties of aqueous sodium chloride solutions,
      ! Journal of Physical and Chemical Reference Data 13 (1) (1984), pp. 1-102.
      vCoef(I)= & !
      &  Fit%vX(1)/T &
      +  Fit%vX(2)   + Fit%vX(3) *P + Fit%vX(4) *P**2 + Fit%vX(5)*P**3                   &
      +  Fit%vX(6)*LOG(T)                                                                &
      + (Fit%vX(7)   + Fit%vX(8) *P + Fit%vX(9) *P**2 + Fit%vX(10)*P**3)*T               &
      + (Fit%vX(11)  + Fit%vX(12)*P + Fit%vX(13)*P**2                  )*T**2            &
      + (Fit%vX(14)  + Fit%vX(15)*P + Fit%vX(16)*P**2 + Fit%vX(17)*P**3)/(T-227.0D0)     &
      + (Fit%vX(18)  + Fit%vX(19)*P + Fit%vX(20)*P**2 + Fit%vX(21)*P**3)/(680.0D0 - T)
      !
    CASE(5) !"PABALAN-1987"
      vCoef(I)= & !
      &   Fit%vX(1)  + Fit%vX(2)*P            &
      + ( Fit%vX(3)  + Fit%vX(4)*P )  /T      &
      +   Fit%vX(5)                   *LOG(T) &
      + ( Fit%vX(6)  + Fit%vX(7)*P )  *T      &
      + ( Fit%vX(8)  + Fit%vX(9)*P )  *T**2   &
      +   Fit%vX(10)                  /(T-227.0D0) &
      + ( Fit%vX(11) + Fit%vX(12)*P ) /(647.0D0-T)
      !
    CASE(6) !"POLYA-2007"
      vCoef(I)= & !
      & Fit%vX(1)    &
      + Fit%vX(2) *T &
      + Fit%vX(3) /T &
      + Fit%vX(4) /(T-210.0D0) &
      + Fit%vX(5) /(647.0D0-T) &
      + Fit%vX(6) *(T-443.0D0)**3/3.0D0 &
      + Fit%vX(7) *(P-One) &
      + Fit%vX(8) *(P-One)**2/2.0D0

    CASE(7) ! "LI-DUAN-2007"
      vCoef(I)= & !
      & Fit%vX(1)             &
      + Fit%vX(2) *T          &
      + Fit%vX(3) /T          &
      + Fit%vX(4)*T**2        &
      + Fit%vX(5)/(630.0D0 - T)   &
      + Fit%vX(6) *P          &
      + Fit%vX(7)*P*LOG(T)    &
      + Fit%vX(8)*P/T         &
      + Fit%vX(9) *P/(630.0D0 - T) &
      + Fit%vX(10)*P**2/(630.0D0 - T)**2  &
      + Fit%vX(11)*T*LOG(P)

    CASE DEFAULT
      CALL Stop_("in Fitt_TPUpdate, iModel out of range")

    END SELECT
    !
  END DO
  !
  RETURN
END SUBROUTINE Fitt_TPUpdate

SUBROUTINE FGENERIQUE (INDIC,NOM,T,P,A,AV,TREF1,TREF2,TREF3,FF)
!!! UNUSED !!!
! Calcul des parametres de Pitzer en fonction de T et P
!
! Arguments d'entree :
! --------------------
! INDIC : si 0 : fonction par defaut (Monnin, Chemical Geology, 153(1999)187-209
!   sinon : l'user specifie sa fonction
! NOM : nom de la fonction si l'user la specifie
! T : Temperature (K)
! P : Pression (Pa)
! A : parametres de la fonction par defaut
! AV: parametres de la fonction par defaut
! TREF1 : temperature de reference 1 pour la fonction par defaut (K)
! TREF2 : temperature de reference 2 pour la fonction par defaut (K)
! TREF3 : temperature de reference 3 pour la fonction par defaut (K)
!
! Arguments de sortie :
! ---------------------
! FF : valeur calculee par la fonction
  !
  INTEGER,      INTENT(IN):: INDIC
  CHARACTER*200,INTENT(IN):: NOM
  REAL(dp),     INTENT(IN):: A(10),AV(10)
  REAL(dp),     INTENT(IN):: T,P,TREF1,TREF2,TREF3
  !
  REAL(dp),     INTENT(OUT):: FF
  !
  REAL(dp):: F0,FV,PBAR
  !
  PBAR=P*1.0d-5
  !
  IF (INDIC==0) THEN
    !----------------------------- Fonction par defaut (Monnin, 1999) --
    F0= A(1)                 &
    & + A(2) *T             &
    & + A(3) /T             &
    & + A(4) *LOG(T)        &
    & + A(5) /(T - TREF1)   &
    & + A(6) *T**2          &
    & + A(7) /(TREF2 - T)   &
    & + A(8) /(T - TREF3)   &
    & + A(9) *T**3          &
    & + A(10)*T**4
    !
    !-- Volume function (pressure correction)
    FV= AV(1)              &
    & + AV(2) *T            &
    & + AV(3) /T            &
    & + AV(4) *LOG(T)       &
    & + AV(5) /(T - TREF1)  &
    & + AV(6) *T**2         &
    & + AV(7) /(TREF2 - T)  &
    & + AV(8) /(T - TREF3)  &
    & + AV(9) *T**3         &
    & + AV(10)*T**4
    !
    FF=  F0 + PBAR * FV
    !----------------------------/ Fonction par defaut (Monnin, 1999) --
  ELSE
    SELECT CASE (NOM)
      !  Ecrire ici CASE ('Nom de la nouvelle fonction')
      !  Ecrire ici FF=
      !  CASE (
    END SELECT
  ENDIF
  !
  RETURN
END SUBROUTINE FGENERIQUE

!~ SUBROUTINE Fitt_TPUpdate_(TdgK,vFitt,vCoef)
  !~ REAL(dp),INTENT(IN):: TdgK
  !~ TYPE(T_Fitt),INTENT(IN):: vFitt(:)
  !~ REAL(dp),INTENT(OUT):: vCoef(:)
  !~ !
  !~ INTEGER:: I
  !~ REAL(dp),PARAMETER:: Tref= 298.15D0
  !~ !
  !~ DO I=1,SIZE(vCoef)
    !~ vCoef(I)= vFitt(I)%X0 &
    !~ &       + vFitt(I)%X1 *(One/TdgK -One/Tref) &
    !~ &       + vFitt(I)%X2 *LOG(TdgK/Tref)       &
    !~ &       + vFitt(I)%X3 *(TdgK     -Tref)     &
    !~ &       + vFitt(I)%X4 *(TdgK*TdgK -Tref*Tref)
  !~ END DO
  !~ !
  !~ RETURN
!~ END SUBROUTINE Fitt_TPUpdate_

SUBROUTINE Solmodel_Pitzer_Dtb_TPtest
  USE M_IoTools
  USE M_Files,ONLY: DirOut
  !
  INTEGER :: f1,f2
  ! INTEGER :: I,J,K
  REAL(dp):: vT(1:9),Pbar
  ! TYPE(T_Fitt):: Fitt
  ! REAL(dp),ALLOCATABLE:: Res(:,:)
  ! CHARACTER(LEN=30):: cFormat
  !
  vT(1:9)=(/ 0.0D0,  25.0D0, 50.0D0, 75.0D0, &
  &        100.0D0, 125.0D0,150.0D0,175.0D0, &
  &        200.0D0 /)
  Pbar= 1.0D0
  !
  !WRITE(cFormat,'(A,I3,A)') '(3(A,1X),',n,'(G12.3,1X))'
  !WRITE(cFormat,'(A)') '(3(A,1X),5(G12.3,1X))'
  !
  !WRITE(ff,cFormat,ADVANCE="no") "X=",(x(k),k=1,n)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Pitzer_Dtb_TPtest"
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirOut)//"_pitzer_tptest.tab")
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirOut)//"_pitzer_tpcoeffs.tab")
  !
  CALL TPtest(f1,f2, "BETA0", vT,Pbar, tI_Beta0, vF_Beta0,vBeta0)
  CALL TPtest(f1,f2, "BETA1", vT,Pbar, tI_Beta1, vF_Beta1,vBeta1)
  CALL TPtest(f1,f2, "BETA2", vT,Pbar, tI_Beta2, vF_Beta2,vBeta2)
  CALL TPtest(f1,f2, "CPHI",  vT,Pbar, tI_CPhi,  vF_CPhi, vCPhi )
  CALL TPtest(f1,f2, "THETA", vT,Pbar, tI_Theta, vF_Theta,vTheta)
  CALL TPtest(f1,f2, "LAMBDA",vT,Pbar, tI_Lamda, vF_Lamda,vLamda)
  !~ ALLOCATE(Res(9,SIZE(vBeta0)))
  !~ DO I=1,9
    !~ CALL Fitt_TPUpdate(vT(I)+273.15D0,Pbar, vF_Beta0, vBeta0)
    !~ Res(I,:)= vBeta0(:)
  !~ ENDDO
  !~ DO I=1,SIZE(tI_Beta0,1)
    !~ DO J=1,SIZE(tI_Beta0,2)

      !~ IF(tI_Beta0(I,J)>0) THEN
        !~ PRINT *, &
        !~ & TRIM(vSolut(I)%NamSp), " ", &
        !~ & TRIM(vSolut(J)%NamSp)
        !~ Fitt= vF_Beta0(tI_Beta0(I,J))
        !~ WRITE(f1,'(3(A12,1X),9(G15.6,1X))') &
        !~ & "BETA0", &
        !~ & TRIM(vSolut(I)%NamSp), &
        !~ & TRIM(vSolut(J)%NamSp), &
        !~ & (Res(K,tI_Beta0(I,J)), K=1,9)
        !~ WRITE(f2,'(4(A12,1X),24(G15.6,1X))') &
        !~ & TRIM(Fitt%namModel), &
        !~ & "BETA0", &
        !~ & TRIM(vSolut(I)%NamSp), &
        !~ & TRIM(vSolut(J)%NamSp), &
        !~ & (Fitt%vX(K),           K=1,Fitt%N)
      !~ ENDIF

    !~ ENDDO
  !~ ENDDO
  !~ DEALLOCATE(Res)
  !
  !~ ALLOCATE(Res(9,SIZE(vBeta1)))
  !~ DO I=1,9
    !~ CALL Fitt_TPUpdate(vT(I)+273.15D0,Pbar, vF_beta1, vbeta1)
    !~ Res(I,:)= vbeta1(:)
  !~ ENDDO
  !~ DO I=1,SIZE(tI_beta1,1)

    !~ DO J=1,SIZE(tI_beta1,2)

      !~ IF(tI_beta1(I,J)>0) THEN
        !~ Fitt= vF_Beta1(tI_Beta1(I,J))
        !~ WRITE(f1,'(3(A12,1X),9(G15.6,1X))') &
        !~ & "BETA1", &
        !~ & TRIM(vSolut(I)%NamSp), &
        !~ & TRIM(vSolut(J)%NamSp), &
        !~ & (Res(K,tI_beta1(I,J)), K=1,9)
        !~ WRITE(f2,'(4(A12,1X),24(G15.6,1X))') &
        !~ & TRIM(Fitt%namModel), &
        !~ & "BETA1", &
        !~ & TRIM(vSolut(I)%NamSp), &
        !~ & TRIM(vSolut(J)%NamSp), &
        !~ & (Fitt%vX(K),           K=1,Fitt%N)
      !~ ENDIF

    !~ ENDDO

  !~ ENDDO
  !~ DEALLOCATE(Res)
  !
  CLOSE(f1)
  CLOSE(f2)
  !
  IF(iDebug>0) PRINT '(/,A,/)',"=!= results in pitzer_tptest.tab =!="
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "<  Solmodel_Pitzer_Pitzer_Dtb_TPtest"
  !
END SUBROUTINE Solmodel_Pitzer_Dtb_TPtest

SUBROUTINE TPtest(f1,f2,Str,vT,Pbar,tI_Beta,vF_Beta,vBeta)
  !
  INTEGER,     INTENT(IN):: f1,f2
  CHARACTER(*),INTENT(IN):: Str
  REAL(dp),    INTENT(IN):: vT(:)
  REAL(dp),    INTENT(IN):: Pbar
  INTEGER,     INTENT(IN):: tI_Beta(:,:)
  TYPE(T_Fitt),INTENT(IN):: vF_Beta(:)
  REAL(dp),    INTENT(IN):: vBeta(:)
  !
  TYPE(T_Fitt):: Fitt
  INTEGER:: I,J,K,nT
  REAL(dp),ALLOCATABLE:: Res(:,:)
  REAL(dp),ALLOCATABLE:: vX(:)
  
  nT= SIZE(vT)
  !
  ALLOCATE(vX(SIZE(vBeta)))
  ALLOCATE(Res(9,SIZE(vBeta)))
  !
  DO I=1,nT
    CALL Fitt_TPUpdate(vT(I)+273.15D0,Pbar, vF_Beta, vX)
    Res(I,:)= vX(:)
  ENDDO
  !
  DO I=1,SIZE(tI_Beta,1)
    DO J=1,SIZE(tI_Beta,2)

      IF(tI_Beta(I,J)>0) THEN
        PRINT *, &
        & TRIM(vSolut(I)%NamSp), " ", &
        & TRIM(vSolut(J)%NamSp)
        Fitt= vF_Beta(tI_Beta(I,J))
        WRITE(f1,'(2(A,1X),9(G15.6,1X))') &
        & TRIM(Str), &
        & TRIM(vSolut(I)%NamSp)//"="//TRIM(vSolut(J)%NamSp), &
        & (Res(K,tI_Beta(I,J)), K=1,nT)
        WRITE(f2,'(3(A,1X),24(G15.6,1X))') &
        & TRIM(Fitt%namModel), &
        & TRIM(Str), &
        & TRIM(vSolut(I)%NamSp)//"="//TRIM(vSolut(J)%NamSp), &
        & (Fitt%vX(K), K=1,Fitt%N)
      ENDIF

    ENDDO
  ENDDO
  !
  DEALLOCATE(Res)
  DEALLOCATE(vX)

END SUBROUTINE TPtest

SUBROUTINE ReadRVals4(Line,x1,x2,x3,x4)
  USE M_IOTools,ONLY:LinToWrd,WrdToReal
  CHARACTER(LEN=*) Line
  CHARACTER(255)   Word
  REAL(dp)::x1,x2,x3,x4
  INTEGER  ::i
  LOGICAL  ::EoL
  REAL(dp),DIMENSION(1:4)::vX
  vX=0.0
  EoL=.TRUE.
  DO i=1,4
    CALL LinToWrd(Line,Word,EoL)
    CALL WrdToReal(Word,vX(i))
    IF(EoL) EXIT
  ENDDO
  x1=vX(1) ; x2=vX(2) ; x3=vX(3) ; x4=vX(4)
ENDSUBROUTINE ReadRVals4

SUBROUTINE Solmodel_Pitzer_Dtb_Check(ff,vSpc)
  USE M_T_Species,ONLY: T_Species
  
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  INTEGER,INTENT(IN):: ff
  
  CALL Dtb_Check(ff,"BETA0",    vSpc, tI_Beta0 , vBeta0)
  CALL Dtb_Check(ff,"BETA1",    vSpc, tI_Beta1 , vBeta1)
  CALL Dtb_Check(ff,"BETA2",    vSpc, tI_Beta2 , vBeta2)
  CALL Dtb_Check(ff,"ALFA1",    vSpc, tI_Beta1 , vAlfa1)
  CALL Dtb_Check(ff,"ALFA2",    vSpc, tI_Beta2 , vAlfa2)
  CALL Dtb_Check(ff,"CPHI",     vSpc, tI_CPhi  , vCPhi)
  CALL Dtb_Check(ff,"THETA",    vSpc, tI_Theta , vTheta)
  CALL Dtb_Check(ff,"LAMBDA",   vSpc, tI_Lamda , vLamda)
  CALL Dtb_Check_Dim3(ff,"ZETA",vSpc, tI_Zeta  , vZeta)
  CALL Dtb_Check_Dim3(ff,"PSI", vSpc, tI_Psi   , vPsi)
  
  RETURN
END SUBROUTINE Solmodel_Pitzer_Dtb_Check

SUBROUTINE Solmodel_Pitzer_Dtb_Check_Old(ff,vSpc)
  USE M_T_Species,ONLY: T_Species
  
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  INTEGER,INTENT(IN):: ff
  
  CALL Dtb_Check_Old(ff,"BETA0",vSpc,tBeta0)
  CALL Dtb_Check_Old(ff,"BETA1",vSpc,tBeta1)
  CALL Dtb_Check_Old(ff,"BETA2",vSpc,tBeta2)
  CALL Dtb_Check_Old(ff,"ALFA1",vSpc,tAlfa1)
  CALL Dtb_Check_Old(ff,"ALFA2",vSpc,tAlfa2)
  CALL Dtb_Check_Old(ff,"CPHI", vSpc,tCPhi)
  CALL Dtb_Check_Old(ff,"THETA",vSpc,tTheta)
  CALL Dtb_Check_Old(ff,"LAMBDA",vSpc,tLamda)
  CALL Dtb_Check_Dim3_Old(ff,"ZETA",vSpc,tZeta)
  CALL Dtb_Check_Dim3_Old(ff,"PSI", vSpc,tPsi)
  
  RETURN
END SUBROUTINE Solmodel_Pitzer_Dtb_Check_Old

SUBROUTINE Dtb_Check(f,Str,vSpc,tI_Coef,vCoef)

  USE M_T_Species,ONLY: T_Species
  !
  INTEGER,        INTENT(IN):: f
  CHARACTER(*),   INTENT(IN):: Str
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  INTEGER,        INTENT(IN):: tI_Coef(:,:)
  REAL(dp),       INTENT(IN):: vCoef(:)
  !
  INTEGER :: I,J
  REAL(dp):: X
  !
  WRITE(f,'(/,A,/)') TRIM(Str)
  DO I=1,SIZE(tI_Coef,1)
    DO J=1,SIZE(tI_Coef,2)
      IF(tI_Coef(I,J)>0) THEN
        X= vCoef(tI_Coef(I,J))
        WRITE(f,'(G15.6,A1,2(A,A1))') &
        & X,t_,TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_
      ENDIF
    ENDDO
  ENDDO
  !
ENDSUBROUTINE Dtb_Check

SUBROUTINE Dtb_Check_Dim3(f,Str,vSpc,tI_Coef,vCoef)
!--
!-- same as Dtb_Check, for 3D tables
!--
  USE M_T_Species,ONLY: T_Species
  !
  INTEGER,        INTENT(IN):: f
  CHARACTER(*),   INTENT(IN):: Str
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  INTEGER,        INTENT(IN):: tI_Coef(:,:,:)
  REAL(dp),       INTENT(IN):: vCoef(:)
  !---------------------------------------------------------------------
  INTEGER :: I,J,K
  REAL(dp):: X
  !
  WRITE(f,'(/,A,/)') TRIM(Str)
  DO I=1,SIZE(tI_Coef,1)
    DO J=1,SIZE(tI_Coef,2)
      DO K=1,SIZE(tI_Coef,3)
        IF(tI_Coef(I,J,K)>0) THEN
          X= vCoef(tI_Coef(I,J,K))
          WRITE(f,'(G15.6,A1,3(A,A1))') &
          & X,t_, &
          & TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_,TRIM(vSpc(K)%NamSp),t_
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  
  RETURN
ENDSUBROUTINE Dtb_Check_Dim3

SUBROUTINE Dtb_Check_Old(f,Str,vSpc,vCoef)
!-- 
  USE M_T_Species,ONLY: T_Species
  !---------------------------------------------------------------------
  INTEGER,        INTENT(IN):: f
  CHARACTER(*),   INTENT(IN):: Str
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  REAL(dp),       INTENT(IN):: vCoef(:,:)
  !
  INTEGER :: I,J
  REAL(dp):: X
  !---------------------------------------------------------------------
  
  WRITE(f,'(/,A,/)') TRIM(Str)
  DO I=1,SIZE(vCoef,1)
    DO J=1,SIZE(vCoef,2)
      X= vCoef(I,J)
      IF(X/=Zero) &
      & WRITE(f,'(G15.6,A1,2(A,A1))') &
      & X,t_,TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_
    ENDDO
  ENDDO
  
  RETURN
ENDSUBROUTINE Dtb_Check_Old

SUBROUTINE Dtb_Check_Dim3_Old(f,Str,vSpc,vCoef)
! same as Dtb_Check, for 3D tables
  USE M_T_Species,ONLY: T_Species
  !
  INTEGER,        INTENT(IN):: f
  CHARACTER(*),   INTENT(IN):: Str
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  REAL(dp),       INTENT(IN):: vCoef(:,:,:)
  !
  INTEGER :: I,J,K
  REAL(dp):: X
  !
  WRITE(f,'(/,A,/)') TRIM(Str)
  DO I=1,SIZE(vCoef,1)
    DO J=1,SIZE(vCoef,2)
      DO K=1,SIZE(vCoef,3)
        X= vCoef(I,J,K)
        IF(X/=Zero) &
        & WRITE(f,'(G15.6,A1,3(A,A1))') &
        & X,t_,TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_,TRIM(vSpc(K)%NamSp),t_
      ENDDO
    ENDDO
  ENDDO
  !
ENDSUBROUTINE Dtb_Check_Dim3_Old

SUBROUTINE Solmodel_Pitzer_Dtb_Clean
  !
  IF (ALLOCATED(vSolut)) DEALLOCATE(vSolut)
  !
  IF (ALLOCATED(tBeta0))   DEALLOCATE(tBeta0)
  IF (ALLOCATED(tBeta1))   DEALLOCATE(tBeta1)
  IF (ALLOCATED(tBeta2))   DEALLOCATE(tBeta2)
  IF (ALLOCATED(tCPhi))    DEALLOCATE(tCPhi)
  IF (ALLOCATED(tTheta))   DEALLOCATE(tTheta)
  IF (ALLOCATED(tLamda))   DEALLOCATE(tLamda)
  IF (ALLOCATED(tAlfa1))   DEALLOCATE(tAlfa1)
  IF (ALLOCATED(tAlfa2))   DEALLOCATE(tAlfa2)
  IF (ALLOCATED(tZeta))    DEALLOCATE(tZeta)
  IF (ALLOCATED(tPsi))     DEALLOCATE(tPsi)
  !
  IF (ALLOCATED(tI_Beta0))   DEALLOCATE(tI_Beta0)
  IF (ALLOCATED(tI_Beta1))   DEALLOCATE(tI_Beta1)
  IF (ALLOCATED(tI_Beta2))   DEALLOCATE(tI_Beta2)
  IF (ALLOCATED(tI_CPhi))    DEALLOCATE(tI_CPhi)
  IF (ALLOCATED(tI_Theta))   DEALLOCATE(tI_Theta)
  IF (ALLOCATED(tI_Lamda))   DEALLOCATE(tI_Lamda)
  IF (ALLOCATED(tI_Zeta))    DEALLOCATE(tI_Zeta)
  IF (ALLOCATED(tI_Psi))     DEALLOCATE(tI_Psi)
  !
  IF (ALLOCATED(vF_Beta0))   DEALLOCATE(vF_Beta0)
  IF (ALLOCATED(vF_Beta1))   DEALLOCATE(vF_Beta1)
  IF (ALLOCATED(vF_Beta2))   DEALLOCATE(vF_Beta2)
  IF (ALLOCATED(vF_CPhi))    DEALLOCATE(vF_CPhi)
  IF (ALLOCATED(vF_Theta))   DEALLOCATE(vF_Theta)
  IF (ALLOCATED(vF_Lamda))   DEALLOCATE(vF_Lamda)
  IF (ALLOCATED(vF_Zeta))    DEALLOCATE(vF_Zeta)
  IF (ALLOCATED(vF_Psi))     DEALLOCATE(vF_Psi)
  !
  IF (ALLOCATED(vBeta0))   DEALLOCATE(vBeta0)
  IF (ALLOCATED(vBeta1))   DEALLOCATE(vBeta1)
  IF (ALLOCATED(vBeta2))   DEALLOCATE(vBeta2)
  IF (ALLOCATED(vAlfa1))   DEALLOCATE(vAlfa1)
  IF (ALLOCATED(vAlfa2))   DEALLOCATE(vAlfa2)
  IF (ALLOCATED(vCPhi))    DEALLOCATE(vCPhi)
  IF (ALLOCATED(vTheta))   DEALLOCATE(vTheta)
  IF (ALLOCATED(vLamda))   DEALLOCATE(vLamda)
  IF (ALLOCATED(vZeta))    DEALLOCATE(vZeta)
  IF (ALLOCATED(vPsi))     DEALLOCATE(vPsi)
  !
ENDSUBROUTINE Solmodel_Pitzer_Dtb_Clean

ENDMODULE M_Solmodel_Pitzer_Dtb
