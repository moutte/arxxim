MODULE M_KinModel_Read
!--
!-- routines for reading kinetic models
!--
  USE M_Kinds
  USE M_Trace,     ONLY: Stop_,iDebug,fTrc,T_, Fatal_, Warning_,Pause_
  USE M_T_Kinmodel,ONLY: T_KinModel
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: KinModel_Init
  !
  !PUBLIC:: KinModel_FileToLnk
  !PUBLIC:: KinModel_Alloc
  !
  TYPE T_LnkKinModel
    TYPE(T_KinModel)   :: Value
    TYPE(T_LnkKinModel),POINTER:: Next
  ENDTYPE T_LnkKinModel
  TYPE(T_LnkKinModel),POINTER:: Lnk
  !
  !TYPE(T_KinModel),DIMENSION(:),ALLOCATABLE:: vKinModel
  !
CONTAINS

SUBROUTINE KinModel_Init(vSpc)
  USE M_T_Species,  ONLY: T_Species
  USE M_Global_Vars,ONLY: vKinModel
  
  TYPE(T_Species),DIMENSION(:),INTENT(IN):: vSpc
  
  TYPE(T_LnkKinModel),POINTER:: Lnk
  INTEGER::N, i
  
  CALL KinModel_FileToLnk(vSpc,N,Lnk) !-> available kinetic DATA on minerals
  !
  IF(N>0) THEN
    IF(ALLOCATED(vKinModel)) DEALLOCATE(vKinModel)
    ALLOCATE(vKinModel(1:N))
    CALL KinModel_Alloc(Lnk,vKinModel)
  ELSE
    IF(iDebug>0) WRITE(fTrc,'(A)') "NO Kinetic Data Found !!!!!!!!!!!!!!"
    IF(iDebug>2) PRINT *,"NO Kinetic Database Found !!!"
  ENDIF
  !
  IF(iDebug==4) THEN
    PRINT *,"KinModel_Init >>"
    DO i=1,SIZE(vKinModel)
      PRINT *,vKinModel(i)%Name
    ENDDO
    CALL pause_
  ENDIF
  
  RETURN
ENDSUBROUTINE KinModel_Init

SUBROUTINE BuildLnk(B,E,Lnk,pCur)
  LOGICAL            :: B !TRUE=first element
  TYPE(T_KinModel)   :: E
  TYPE(T_LnkKinModel),POINTER:: Lnk,pCur
  
  IF(B) NULLIFY(Lnk)
  
  IF(B) THEN
    ALLOCATE(Lnk)
    NULLIFY(Lnk%next)
    Lnk%Value=E
    pCur => Lnk
  ELSE
    ALLOCATE(pCur%next)
    NULLIFY(pCur%next%next)
    pCur%next%Value=E
    pCur => pCur%next
  ENDIF
  
ENDSUBROUTINE BuildLnk

SUBROUTINE KinModel_FileToLnk( &
& vSpc, & !IN
& N,Lnk)  !OUT
!--
!-- read kinetic data for minerals (and any other species)
!-- -> build LnkKin
!--
  USE M_IOTools !, ONLY:dimV,LinToWrd,ReadRValsV
  USE M_Files,   ONLY: NamFKin,DirLog
  USE M_Global_Vars,  ONLY: nAq !,vSpc,nSp
  USE M_T_Species, ONLY: T_Species,Species_Index,Species_Rename
  USE M_T_Kinmodel,ONLY: MaxKinTerm,T_KinModel
  !
  TYPE(T_Species),DIMENSION(:),INTENT(IN):: vSpc
  INTEGER,            INTENT(OUT):: N
  TYPE(T_LnkKinModel),POINTER    :: Lnk
  !
  TYPE(T_LnkKinModel),POINTER::p
  CHARACTER(LEN=512):: L,W,W1
  TYPE(T_KinModel)  ::MK
  LOGICAL :: EoL,NewM,OldFormat
  LOGICAL :: LPrecip,LDissol,ModelIsOk
  INTEGER :: iAq,J
  INTEGER :: f,ios
  REAL(dp):: X1,X2,X3
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< KinModel_FileToLnk"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFKin))
  !
  OldFormat=  .FALSE.
  N= 0
  !
  DoFile: DO 
    
    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    
    S1: SELECT CASE(W)
    
    CASE("ENDINPUT") S1
      EXIT  DoFile  
    
    CASE("KINETICS") S1
    !format for kinetic models, a la USGS, new version with possibly two species
      !
      IF(.NOT. EoL) THEN
        CALL LinToWrd(L,W1,EoL)
        OldFormat= TRIM(W1)=="OLD"
      ENDIF
      !
      IF(OldFormat) THEN
      !
      N=0; NewM=.FALSE.
      !
      DoBlockOld: DO !build LnkKin of all "kinetic" minerals
      
        READ(F,'(A)',IOSTAT=ios) L
        IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoBlockOld !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        
        S2_: SELECT CASE(W)
          CASE("ENDINPUT") S2_; EXIT DoFile !to prevent bad file format 
          CASE("END","ENDKINETICS") S2_; EXIT DoBlockOld
        END SELECT S2_
        
        IF(W(1:1)/='&') THEN
          !-- IF first char is not '&', the line contain a name
          IF(NewM) THEN
            CALL BuildLnk(N==1,MK,Lnk,p)
            NewM=.FALSE.
            !!IF(iDebug>0) WRITE(fTrc,'(A24,I3,A24)') "Min=',M%Name,I," =Spc ",Species_Index(M%Name)
          ENDIF
          !
          !!MK%Special=0
          MK%NTermD=0; MK%pK_d=Zero;  MK%E_d=Zero; MK%N_d=Zero; MK%AlfaD=One; MK%BetaD=One
          MK%NTermP=0; MK%pK_p=Zero;  MK%E_p=Zero; MK%N_p=Zero; MK%AlfaP=One; MK%BetaP=One
          MK%ISpcDiss(:)=0; MK%JSpcDiss(:)=0 
          MK%ISpcPrec(:)=0; MK%JSpcPrec(:)=0 
          !
          MK%Name=TRIM(W)
          IF(.NOT. EoL) CALL LinToWrd(L,W,EoL)
          J=Species_Index(TRIM(W),vSpc)
          !check whether the mineral is also in the thermodynamic dtb
          IF(J>0) THEN !IF mineral is Ok, THEN set NewM true -> will be saved in LnkKin
            N=N+1
            NewM=.TRUE.
          ENDIF
          
        ELSE
          !-- W(1:1)--'&'
          !-- IF first char is '&',
          !-- the line is the continuation of preceding, i.e. contains data
          IF(NewM) THEN
            !-- READ first word after &, should be either DISSOL or PRECIP
            CALL LinToWrd(L,W,Eol)
            !
            SELECT CASE(TRIM(W))
              CASE("PRECIP")
                LPrecip=.TRUE.; LDissol=.FALSE.
              CASE("DISSOL")
                LDissol=.TRUE.; LPrecip=.FALSE.
              CASE DEFAULT
                CALL Stop_(TRIM(W)//"should have either PRECIP or DISSOL !!!...")
            END SELECT
            !
            CALL LinToWrd(L,W,Eol) !-> W will contain species name, or "QSK"
            !
            IF(TRIM(W)=="QSK") THEN
              !--- read Alpha, Beta
              CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X1)
              IF(LPrecip) MK%AlfaP=X1
              IF(LDissol) MK%AlfaD=X1 
              
              CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X1)
              IF(LPrecip) MK%BetaP=X1
              IF(LDissol) MK%BetaD=X1 
            
            ELSE
              !--- read activator parameters
              CALL Species_Rename(W)
              iAq= Species_Index(TRIM(W),vSpc)
              IF (iAq==0) THEN
                IF(iDebug>0) &
                & WRITE(fTrc,'(2A)') &
                & "in KinModel "//TRIM(MK%Name),", unknown species : "//TRIM(W)
              ENDIF
              !
              CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X1) !-> read kinetic parameters
              CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X2) !-> read kinetic parameters
              CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X3) !-> read kinetic parameters
              !
              IF(iAq>0 .AND. iAq<=nAq) THEN
                !! IF(iDebug>0) WRITE(fTrc,'(A,I3,2A)') "iAq=",iAq," <- ",TRIM(W)
                IF(LPrecip .AND. MK%NTermP<MaxKinTerm) THEN
                !activation energies in kiloJoule in input file !!!!!!!!!!!!!
                  MK%NTermP= MK%NTermP+1
                  MK%ISpcPrec(MK%NTermP)= iAq !-> variable, updated in KinModel_UpDate
                  MK%pK_P(MK%NTermP)=     X1
                  MK%E_P(MK%NTermP)=      X2*1.D3
                  MK%N_P(MK%NTermP)=      X3
                ENDIF
                IF(LDissol .AND. MK%NTermD<MaxKinTerm) THEN
                  MK%NTermD= MK%NTermD+1
                  MK%ISpcDiss(MK%NTermD)= iAq !-> variable, updated in KinModel_UpDate
                  MK%pK_D(MK%NTermD)=     X1
                  MK%E_D(MK%NTermD)=      X2*1.D3
                  MK%N_D(MK%NTermD)=      X3
                ENDIF
              ENDIF
              
            ENDIF
            !! M%KModel=MK
          ENDIF !IF NEWM
        ENDIF
      ENDDO DoBlockOld
      
      IF(NewM) CALL BuildLnk(N==1,MK,Lnk,p) !SAVE last mineral in list
      !
      ELSE
      !
      N=0; NewM=.FALSE.
      !
      DoBlock: DO !build LnkKin of all "kinetic" minerals
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        S2: SELECT CASE(W)
          CASE("ENDINPUT") S2; EXIT DoFile !to prevent bad file format 
          CASE("END","ENDKINETICS") S2; EXIT DoBlock
        END SELECT S2
        !
        IF(W(1:1)/='&') THEN
          !-- if first char is not '&', the line contains the model name
          !
          !-- first save model from current buffer to list
          IF(NewM .AND. ModelIsOk) THEN
            N=N+1
            CALL BuildLnk(N==1,MK,Lnk,p)
            NewM=.FALSE.
            !!IF(iDebug>0) WRITE(fTrc,'(A24,I3,A24)') "Min=',M%Name,I," =Spc ",Species_Index(M%Name)
          ENDIF
          !
          MK%Name=TRIM(W)
          !
          MK%NTermD=0; MK%pK_d=Zero;  MK%E_d=Zero; MK%N_d=Zero; MK%AlfaD=One; MK%BetaD=One
          MK%NTermP=0; MK%pK_p=Zero;  MK%E_p=Zero; MK%N_p=Zero; MK%AlfaP=One; MK%BetaP=One
          MK%ISpcDiss(:)=0; MK%JSpcDiss(:)=0 
          MK%ISpcPrec(:)=0; MK%JSpcPrec(:)=0 
          !
          NewM=.TRUE.
          ModelIsOk= .TRUE.
          !
          !QUARTZ USGS2004-1068
          !&  DISSOL  1     1     QsK      
          !&  DISSOL  13.4  90.9  H2O  0    
          !&  PRECIP  1     1     QsK      
          !&  PRECIP  13.4  90.9  H2O  0    
        ELSE
          !-- W(1:1)--'&'
          !-- IF first char is '&',
          !-- the line is the continuation of preceding, i.e. contains data
          IF(NewM) THEN
            CALL LinToWrd(L,W,Eol) !READ first word after &, should be either DISSOL or PRECIP
            SELECT CASE(TRIM(W))
              CASE("PRECIP"); LPrecip=.TRUE.; LDissol=.FALSE.
              CASE("DISSOL"); LDissol=.TRUE.; LPrecip=.FALSE.
              CASE DEFAULT
                CALL Stop_(TRIM(W)//"should have either PRECIP or DISSOL !!!...")
            END SELECT
            !
            CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X1)
            CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X2)
            !
            CALL LinToWrd(L,W,Eol) !-> W will contain a species name, or "QSK"
            !
            IF(TRIM(W)=="QSK") THEN
              !-- read Alpha, Beta
              IF(LPrecip) MK%AlfaP=X1
              IF(LDissol) MK%AlfaD=X1
              IF(LPrecip) MK%BetaP=X2
              IF(LDissol) MK%BetaD=X2
              !
            ELSE
              !-- read activator parameters
              CALL Species_Rename(W)
              iAq=Species_Index(TRIM(W),vSpc)
              !
              IF (iAq==0) THEN
                IF(iDebug>0) &
                & WRITE(fTrc,'(3A)') &
                & "in KinModel "//TRIM(MK%Name), &
                & ", unknown species : "//TRIM(W), &
                & ", >> KINETIC MODEL NOT INCLUDED"
                ModelIsOk= .FALSE.
              ENDIF
              !
              CALL LinToWrd(L,W,Eol)  ;  CALL WrdToReal(W,X3)
              !
              IF(iAq>0 .AND. iAq<=nAq) THEN
                IF(iDebug>0) WRITE(fTrc,'(A,I3,2A)') "iAq=",iAq," <- ",vSpc(iAq)%NamSp
                !
                IF(LPrecip .AND. MK%NTermP<MaxKinTerm) THEN
                  !
                  MK%NTermP= MK%NTermP +1
                  MK%ISpcPrec(MK%NTermP)=iAq
                  !
                  MK%pK_P(MK%NTermP)= X1
                  MK%E_P(MK%NTermP)=  X2*1.D3 !activation energies in kiloJoule in input file
                  MK%N_P(MK%NTermP)=  X3
                  !
                  MK%JSpcPrec(MK%NTermP)=0
                  !
                  IF(.NOT.EoL) THEN !CASE of second species
                    CALL LinToWrd(L,W,Eol)
                    CALL Species_Rename(W)
                    iAq=Species_Index(TRIM(W),vSpc)
                    IF(iAq>0 .AND. iAq<=nAq) THEN 
                      MK%JSpcPrec(MK%NTermP)=iAq
                      CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X1)
                      MK%NJ_P(MK%NTermP)=X1
                    ELSE
                      IF(iDebug>0) &
                      & WRITE(fTrc,'(3A)') &
                      & "in KinModel "//TRIM(MK%Name), &
                      & ", unknown species : "//TRIM(W), &
                      & ", >> KINETIC MODEL NOT INCLUDED"
                      ModelIsOk= .FALSE.
                    ENDIF
                  ENDIF
                  !
                ENDIF
                !
                IF(LDissol .AND. MK%NTermD<MaxKinTerm) THEN
                  !
                  MK%NTermD= MK%NTermD +1
                  MK%ISpcDiss(MK%NTermD)=iAq
                  !
                  MK%pK_D(MK%NTermD)= X1
                  MK%E_D(MK%NTermD)=  X2*1.D3
                  !!! activation energies in kiloJoule in input file !!!
                  MK%N_D(MK%NTermD)=  X3
                  !
                  MK%JSpcDiss(MK%NTermD)=0
                  !
                  IF(.NOT.EoL) THEN !CASE of second species
                    CALL LinToWrd(L,W,Eol)
                    iAq=Species_Index(TRIM(W),vSpc)
                    IF(iAq>0 .AND. iAq<=nAq) THEN 
                      MK%JSpcDiss(MK%NTermD)=iAq
                      CALL LinToWrd(L,W,Eol); CALL WrdToReal(W,X1)
                      MK%NJ_D(MK%NTermD)=X1
                    ENDIF
                  ENDIF
                ENDIF
                !
              ENDIF
              
            ENDIF
          ENDIF !IF NEWM
        ENDIF !IF W(1:1)=='&'
      ENDDO DoBlock
      IF(NewM .AND. ModelIsOk) THEN
        N= N+1
        CALL BuildLnk(N==1,MK,Lnk,p) !save last mineral in list
      ENDIF
      !
      ENDIF
    !ENDCASE("KINETICS")  
    END SELECT S1
     
  ENDDO DoFile
  !
  CLOSE(f)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,I3)') "N=", N
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ KinModel_FileToLnk"
ENDSUBROUTINE KinModel_FileToLnk

SUBROUTINE KinModel_Alloc(Lnk,vKinModel)
  !
  TYPE(T_LnkKinModel),POINTER:: Lnk
  TYPE(T_KinModel),DIMENSION(:),INTENT(OUT):: vKinModel
  !
  TYPE(T_LnkKinModel),POINTER::pCur, pPrev
  INTEGER::I
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< KinModel_Alloc"
  I=0
  pCur=>Lnk
  DO WHILE (ASSOCIATED(pCur))
    I= I+1
    vKinModel(I)=pCur%Value 
    IF(iDebug>0) WRITE(fTrc,'(I4,A1,A12)') I," ",vKinModel(I)%Name
    pPrev=>pCur; pCur=> pCur%next; DEALLOCATE(pPrev)
  ENDDO
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ KinModel_Alloc"
ENDSUBROUTINE KinModel_Alloc

ENDMODULE M_KinModel_Read


