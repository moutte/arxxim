MODULE M_KinFas_Read
!--
!-- routines for reading kinetic models
!--
  USE M_Kinds
  USE M_Trace,   ONLY: iDebug,fTrc,T_,Stop_,Warning_
  USE M_T_KinFas,ONLY: T_KinFas
  USE M_T_Phase, ONLY: T_Phase
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: KinFas_BuildLnk
  PUBLIC:: KinFas_LnkToVec
  PUBLIC:: T_LnkKin 
  !
  TYPE T_LnkKin
    TYPE(T_KinFas)         ::Value
    TYPE(T_LnkKin), POINTER::Next
  ENDTYPE T_LnkKin
  !
CONTAINS

SUBROUTINE KinFas_LnkToVec(LnkKin,vFas,vKinFas)
  !
  TYPE(T_LnkKin),POINTER    :: LnkKin
  TYPE(T_Phase), INTENT(IN) :: vFas(:)
  TYPE(T_KinFas),INTENT(OUT):: vKinFas(:)
  !
  TYPE(T_LnkKin),POINTER::pCur, pPrev
  INTEGER:: I, J
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< KinFas_LnkToVec"
  !
  I=0
  J=0
  pCur=> LnkKin
  DO WHILE (ASSOCIATED(pCur))
    !
    I= I+1
    vKinFas(I)= pCur%Value
    !
    IF(iDebug>0) WRITE(fTrc,'(I4,A1,A12)') I," ",vKinFas(I)%NamKF
    !
    pPrev=>pCur; pCur=> pCur%next; DEALLOCATE(pPrev)
    !
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ KinFas_LnkToVec"
  !
ENDSUBROUTINE KinFas_LnkToVec

SUBROUTINE BuildLnkKin(B,E,L,P)
  LOGICAL               ::B !TRUE=first element
  TYPE(T_KinFas)        ::E
  TYPE(T_LnkKin),POINTER::L,P
  IF(B) NULLIFY(L)
  IF(B) THEN; ALLOCATE(L);      NULLIFY(L%next);      L%Value=E;       P=>L
  ELSE;       ALLOCATE(P%next); NULLIFY(P%next%next); P%next%Value=E;  P=>P%next
  ENDIF
ENDSUBROUTINE BuildLnkKin

SUBROUTINE KinFas_BuildLnk(vFas,vKinModel,sModelSurf,N,LnkKin)
!--
!-- read the "ROCK" block from input.file -> build link list --
!--
  USE M_IOTools
  USE M_Files,      ONLY: NamFInn
  !!USE M_T_Species, ONLY: T_Species,Species_Index
  USE M_T_Phase,    ONLY: T_Phase,Phase_Index
  USE M_T_Kinmodel, ONLY: T_KinModel,KinModel_Index
  USE M_T_KinFas,   ONLY: KinFas_Zero,KinFas_Surf_Zero
  !
  TYPE(T_Phase),   INTENT(IN) :: vFas(:)
  TYPE(T_KinModel),INTENT(IN) :: vKinModel(:)
  CHARACTER(LEN=7),INTENT(IN) :: sModelSurf
  INTEGER,         INTENT(OUT):: N
  TYPE(T_LnkKin),  POINTER    :: LnkKin
  !
  REAL(dp),DIMENSION(dimV)::vX
  TYPE(T_LnkKin),POINTER::pKin
  TYPE(T_KinFas)    :: M
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W1,W2,W3
  REAL(dp):: xGram,xMole
  LOGICAL :: EoL,BlockFound,OldFormat
  INTEGER :: I,jKin,mDum,f_,ios
  !
  !dimV: defined in ModIo,
  !dimension of "buffer array" scanned from words of a string
  !
  CALL GetUnit(f_)
  OPEN(f_,FILE=TRIM(NamFInn))
  !
  N= 0
  !
  BlockFound= .FALSE.
  OldFormat=  .FALSE.
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< KinFas_BuildLnk"
  DoFile: DO 
    !
    READ(F_,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W1,EoL)
    IF(W1(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W1,EoL)
    !
    SELECT CASE(W1)
    !
    CASE("ENDINPUT"); EXIT DoFile
    !
    !NB: the DYNAMIC.ROCK block is READ after KINETICS block
    !
    CASE("DYNAMIC.ROCK")
      BlockFound=.TRUE.
      N=0
      IF(.NOT. EoL) THEN
        CALL LinToWrd(L,W3,EoL)
        OldFormat= TRIM(W3)=="OLD"
      ENDIF
      IF(OldFormat) THEN
      !------------------------------------------------------ old format
        DoReadOld: DO
        
          READ(F_,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W1,EoL)
          IF(W1(1:1)=='!')   CYCLE DoReadOld !skip comment lines
          CALL AppendToEnd(L,W1,EoL)
          !
          SELECT CASE(W1)
          CASE("ENDINPUT")                ; EXIT DoFile
          CASE("END","ENDDYNAMIC.ROCK")   ; EXIT DoReadOld
          END SELECT
          !
          CALL LinToWrd(L,W2,EoL)
          !
          IF(TRIM(W2)/="D" .AND. TRIM(W2)/="R") THEN
            !-> IF next word is not D or R, THEN it is Species Name
            
            !I= Species_Index(W2,vSpc) !find species named W2 in vSpc
            I= Phase_Index(W2,vFas) !find phase named W2 in vFas
            
            CALL LinToWrd(L,W2,EoL) !-> next word is Name Of Kinetic Model
            jKin= KinModel_Index(W2,vKinModel) !find species named W2 in vKinModel
            
            CALL LinToWrd(L,W2,EoL) !READ next word
            IF(TRIM(W2)/="D" .AND. TRIM(W2)/="R") CALL Stop_( &
            & "texture code should be either R (radius) or D (density)" )
            M%cMode=W2(1:1)
          
          ELSE
            I= Phase_Index(W1,vFas)
            !I= Species_Index(W1,vSpc)
            jKin= KinModel_Index(W1,vKinModel)
            M%cMode=W2(1:1)
          
          ENDIF
          !
          IF(I==0) WRITE(*,'(A)') "Species not found for "//TRIM(W1)
          IF(jKin==0) WRITE(*,'(A)') "Kinetic model not found for "//TRIM(W1)
          !
          IF(I>0 .AND. jKin>0) THEN
            !
            CALL KinFas_Zero(M)
            !
            M%NamKF=TRIM(W1) !-> the Kinetic Mineral Name
            !M%iSpc=I
            M%iFas=I
            M%iKin=jKin
            !
            IF(iDebug>0) WRITE(fTrc,'(3(A12,1X))') &
            !& M%Name,vSpc(M%iSpc)%Name,vKinModel(M%iKin)%Name
            & M%NamKF,vFas(M%iFas)%NamFs,vKinModel(M%iKin)%Name
            !
            CALL ReadRValsV(L,mDum,vX) !read radius and relative fraction
            M%Dat%Radius= vX(1)
            M%Dat%PhiM=   vX(2)
            IF(vX(3)/=Zero) M%QsKSeuil= vX(3)
            IF(vX(4)/=Zero) M%ReacCoef= vX(4)
            !
            !!M%Dat%cSat="1" !"PRIMARY" 
            !!IF(M%Dat%PhiM<1.0E-5) M%Dat%cSat="2" !"SECONDA" 
            !!!-> mineral is "secondary",i.e. not REALly in initial assemblage
            !
            N=N+1; CALL BuildLnkKin(N==1,M,LnkKin,pKin)
          ENDIF
        ENDDO DoReadOld
        !
      !------------------------------------------------------/old format
      ELSE
      !------------------------------------------------------ new format
        !
        DoReadRock: DO
        
          READ(F_,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W1,EoL)
          IF(W1(1:1)=='!')   CYCLE DoReadRock !skip comment lines
          CALL AppendToEnd(L,W1,EoL)
          SELECT CASE(W1)
            CASE("ENDINPUT")              ; EXIT DoFile
            CASE("END","ENDDYNAMIC.ROCK") ; EXIT DoReadRock
          END SELECT
          !
          I=    0
          jKin= 0
          xGram= Zero
          xMole= Zero
          !
          CALL KinFas_Zero(M)
          !
          M%NamKF= TRIM(W1) !-> the Kinetic Phase Name
          CALL LinToWrd(L,W2,EoL)
          !
          DO
            
            IF(EoL) EXIT
            !
            SELECT CASE(TRIM(W2))
            !
            CASE DEFAULT
              CALL Stop_(TRIM(W2)//" <<Unknown Keyword in DYNAMIC.ROCK")
            !
            CASE("SPECIES","SOLUTION","MIXTURE")
              IF(TRIM(W2)=="SOLUTION") &
              & CALL Warning_(TRIM(W2)//" soon Obsolete, better use MIXTURE !!!")
              !
              CALL LinToWrd(L,W2,EoL)
              I= Phase_Index(W2,vFas)
              IF(I==0) PRINT '(A)',"SPECIES not found for "//TRIM(W2)
            !
            CASE("KINETICS")
              CALL LinToWrd(L,W2,EoL)
              jKin= KinModel_Index(W2,vKinModel)
              IF(jKin==0) &
              & CALL Stop_("DYNAMIC.ROCK: kinetic model not found for "//TRIM(W2))
            !
            CASE("MODE") !=> mode= Radius / Density / Surface/ Equil / Adsorption
              CALL LinToWrd(L,W2,EoL) ; M%cMode= W2(1:1)
            !
            CASE("RADIUS") !=> radius (metre)
              CALL LinToWrd(L,W2,EoL) ; CALL WrdToReal(W2,M%Dat%Radius)
            !
            CASE("SURFACE") !=> specific surface m2/kg
              CALL LinToWrd(L,W2,EoL) ; CALL WrdToReal(W2,M%Dat%SurfKg) 
            !
            CASE("VOLUME") !=> volume (relative to other volumes in min'assemblages)
              CALL LinToWrd(L,W2,EoL) ; CALL WrdToReal(W2,M%Dat%PhiM)
            !
            CASE("MOLE") !=> mole number
              CALL LinToWrd(L,W2,EoL) ; CALL WrdToReal(W2,xMole)
            !
            CASE("GRAM") !=> gram number
              CALL LinToWrd(L,W2,EoL) ; CALL WrdToReal(W2,xGram)
            !
            END SELECT
            !
            CALL LinToWrd(L,W2,EoL)
            
          ENDDO
          !
          IF(M%Dat%SurfKg >Zero) M%Dat%Radius= -One
          ! otherwise Radius would take default value set in KinFas_Zero ...
          !
          IF(I==0) I= Phase_Index(W1,vFas)
          !-> default species name= name of kinetic phase
          IF(I==0) &
          & CALL Stop_("DYNAMIC.ROCK: Species not found for "//TRIM(M%NamKF))
          !
          IF(M%cMode/="E" .AND. M%cMode/="I") THEN
            IF(jKin==0) jKin= KinModel_Index(W1,vKinModel) !default kin'model name= name of kinetic phase
            IF(jKin==0) CALL Stop_("Kinetic model not found for "//TRIM(M%NamKF))
          ENDIF
          !
          M%iFas= I
          IF(xGram>Zero) M%Dat%PhiM= xGram /1.0D3 /vFas(M%iFas)%WeitKg *vFas(M%iFas)%VolM3
          IF(xMole>Zero) M%Dat%PhiM= xMole *vFas(M%iFas)%VolM3
          !
          IF(M%cMode/="E".AND. M%cMode/="I") THEN
            M%iKin= jKin
          ELSE
            M%iKin= 0
          ENDIF
          !
          CALL KinFas_Surf_Zero(vFas,M)
          !-> calc' SurfKg from Radius, or Radius from SurfKg
          !
          N=N+1
          CALL BuildLnkKin(N==1,M,LnkKin,pKin)
          
        ENDDO DoReadRock
        !
      ENDIF
      !------------------------------------------------------/new format
    !ENDCASE("DYNAMIC.ROCK")
    !
    END SELECT
    !
  ENDDO DoFile
  CLOSE(f_)
  !
  IF(.NOT. BlockFound) PRINT '(A)',"!!!WARNING!!! Block DYNAMIC.ROCK Not Found"
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ KinFas_BuildLnk"
  !
ENDSUBROUTINE KinFas_BuildLnk

SUBROUTINE KinFas_Check(vFas_,vKinFas_)
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_T_Phase,  ONLY: T_Phase,Phase_Index
  !
  TYPE(T_Phase), INTENT(IN):: vFas_(:)
  TYPE(T_KinFas),INTENT(IN):: vKinFas_(:)
  !
  INTEGER ::I,J
  REAL(dp)::RhoM
  
  DO I=1,SIZE(vKinFas_)
    !!J=vKinFas_(I)%iSpc !J= index (in vSpc_) of species named vKinFas(I)%Name
    !!RhoM=vSpc_(J)%WeitKg / vSpc_(J)%V0
    J=vKinFas_(I)%iFas !J= index (in vFas_) of species named vKinFas(I)%Name
    RhoM= vFas_(J)%WeitKg / vFas_(J)%VolM3
    WRITE(fTrc,'(A,I3,A1,A,A15,A1,3(A7,F12.6,A1))') &
    & "iMk=",  I,                                  T_, &
    & "vKinFas=", vFas_(J)%NamFs,                  T_, &
    & "RhoM=", RhoM,                               T_, &
    & "VMol=", vFas_(vKinFas_(I)%iFas)%VolM3,      T_, &
    & "logK=", - vFas_(vKinFas_(I)%iFas)%Grt/Ln10, T_
  ENDDO
  
ENDSUBROUTINE KinFas_Check

ENDMODULE M_KinFas_Read


