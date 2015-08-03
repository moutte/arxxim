MODULE M_DiscretModel_Read
!--
!-- routine for reading block(s) MIXTURE.DISCRETIZE (was SOLUTION.DISCRETIZE)
!-- and build pure species from discretization of mixture phase(s)
!--
  USE M_Kinds
  USE M_Trace,         ONLY: iDebug,fTrc,T_,Stop_,Pause_,Warning_
  USE M_T_DiscretModel,ONLY: T_DiscretModel
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DiscretModel_Read
  !
  !~ TYPE:: T_LnkDiscretModel
    !~ TYPE(T_DiscretModel):: Value
    !~ TYPE(T_LnkDiscretModel),POINTER::Next
  !~ ENDTYPE T_LnkDiscretModel
  !
CONTAINS

SUBROUTINE DiscretModel_Read(vMixModel)
!--
!-- called by DiscretPhase_Add in MODULE m_global_alloc
!--
  USE M_T_MixModel,ONLY: T_MixModel
  !
  USE M_Global_Vars,ONLY: vDiscretModel
  !
  TYPE(T_MixModel),INTENT(IN) :: vMixModel(:)  
  !
  TYPE(T_DiscretModel),ALLOCATABLE:: vTmp(:)
  INTEGER:: N
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_Read"
  !
  ALLOCATE(vTmp(100))
  CALL DiscretModel_ReadFile(vMixModel,vTmp,N)
  !
  IF(N>0) THEN
    DEALLOCATE(vDiscretModel)
    ALLOCATE(vDiscretModel(N))
    vDiscretModel(1:N)= vTmp(1:N)
  ENDIF
  !
  DEALLOCATE(vTmp)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretModel_Read"
  !! IF(iDebug>0) PRINT *,"DiscretModel_Read, N=", N
  !
ENDSUBROUTINE DiscretModel_Read

SUBROUTINE DiscretModel_ReadFile( & !
!--
!-- reads the parameters for the "discretization" of a mixture
!--
& vMixModel,     & !in,  database of solution models
& vDiscretModel, & !out
& N)               !OUT
  USE M_IOTools
  USE M_Files,     ONLY: NamFInn
  USE M_T_MixModel,ONLY: T_MixModel
  USE M_T_DiscretModel
  !
  TYPE(T_MixModel),    INTENT(IN) :: vMixModel(:)
  TYPE(T_DiscretModel),INTENT(OUT):: vDiscretModel(:)
  INTEGER,             INTENT(OUT):: N
  !
  TYPE(T_DiscretModel):: DsModl
  TYPE(T_MixModel)    :: MxModl
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios
  INTEGER:: K,iMix,DD,i,DimMix
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_ReadFile"
  !
  N= 0
  !
  IF(SIZE(vDiscretModel)<1) RETURN
  !
  iMix= 0
  CALL GetUnit(F)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    !
    CALL AppendToEnd(L,W,EoL)
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT DoFile
    !
    CASE("SOLUTION.DISCRETIZE","MIXTURE.DISCRETIZE") !for reading one "block"
      !! MIXTURE.DISCRETIZE
      !!   NAME   BIOT
      !!   NUMBER 9
      !!   MODEL  BIOTITE_MG_IDEAL
      !! ENDSOLUTION.DISCRETIZE
      DoBlock: DO
      
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        
        SELECT CASE(W)
          !
          CASE("ENDINPUT"); EXIT DoFile
          !
          CASE DEFAULT
            CALL Warning_("Unknown (or Obsolete) Keyword: "//TRIM(W))
          !
          CASE("END","ENDSOLUTION.DISCRETIZE","ENDMIXTURE.DISCRETIZE")
            !
            IF(iMix>0) THEN
              DimMix= vMixModel(DsModl%iMix)%NPole
              !
              IF(DimMix>3) THEN
                PRINT *,"only binary or ternary models !!!"
              ELSE
                N= N+1
                !
                IF(DimMix==2) DsModl%Dim2= 0
                IF(DimMix==3) DsModl%Dim2= DsModl%Dim1
                CALL DiscretModel_Dim(DsModl%Dim1,DsModl%Dim2,DsModl%DimTot)
                !
                vDiscretModel(N)= DsModl
                !
                IF(iDebug>2) THEN
                PRINT '(3(A,I4),4A)', &
                & "  Dim1=  ", DsModl%Dim1,   &
                & "  Dim2=  ", DsModl%Dim2,   &
                & "  DimTot=", DsModl%DimTot, &
                & "  MxModl=", TRIM(vMixModel(DsModl%iMix)%Name), &
                & "  DsModl=", TRIM(DsModl%Name)
                ENDIF
                !
                IF(iDebug>0) WRITE(fTrc,'(3(A,I4),4A)') &
                & "  Dim1=  ", DsModl%Dim1,   &
                & "  Dim2=  ", DsModl%Dim2,   &
                & "  DimTot=", DsModl%DimTot, &
                & "  MxModl=", TRIM(vMixModel(DsModl%iMix)%Name), &
                & "  DsModl=", TRIM(DsModl%Name)
              ENDIF
              !
            ENDIF
            EXIT DoBlock
            !
          CASE("NAME")
            CALL LinToWrd(L,W,EoL); DsModl%Name=TRIM(W)
            !
          CASE("NUMBER")
            CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,DD) 
            DD= MIN(DD,MaxDiscret)
            DD= MAX(DD,MinDiscret)
            !
            DsModl%Dim1=   DD -1
            DsModl%DimTot= DsModl%Dim1
            DsModl%Dim2=   0
            !
            IF(.not. Eol) THEN
              CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,DD)
              DD= MIN(DD,MaxDiscret)
              DD= MAX(DD,MinDiscret)
              DsModl%Dim2= DD -1
              DsModl%Dim2= DsModl%Dim1
            ENDIF
            !
          CASE("MODEL")
            CALL LinToWrd(L,W,EoL)
            ! find the solution name in vMixModel%Name
            iMix=0
            DO K=1,SIZE(vMixModel) 
              IF(TRIM(W)==TRIM(vMixModel(K)%Name)) iMix=K
            ENDDO
            !
            IF(iMix>0) THEN
              IF(vMixModel(iMix)%NPole<2) iMix=0
            ENDIF
            IF(iMix==0) CALL Stop_ &
            & ("MIXTURE.DISCRETIZE: "//TRIM(W)//"-> Mixing model not found ...")
            !
            DsModl%iMix= iMix
            DsModl%P1= 1 !default values
            DsModl%P2= 2 !id
            DsModl%P3= 0 !id
            !
          CASE("POLE")
            IF(iMix>0) THEN
              MxModl= vMixModel(iMix)
              DO
                CALL LinToWrd(L,W,EoL)
                I=0
                DO K=1,MxModl%NPole
                  IF(TRIM(W)==TRIM(MxModl%vNamPole(K))) I=K
                ENDDO
                IF(EoL) EXIT
              ENDDO
            ENDIF
            !
        ENDSELECT
        
      ENDDO DoBlock
      
    ENDSELECT
    
  ENDDO DoFile
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretModel_ReadFile"
  !
ENDSUBROUTINE DiscretModel_ReadFile

SUBROUTINE DiscretModel_Dim(Dim1,Dim2,N)
  INTEGER,INTENT(IN) :: Dim1,Dim2
  INTEGER,INTENT(OUT):: N
  !
  INTEGER:: Dims,i0,i,j,k
  !
  Dims= Dim1 +1
  IF(Dim2>0) THEN  ;  i0= 2  ! ternary
  ELSE             ;  i0= 1  ! binary
  ENDIF
  N=0
  i=0
  DO
    i= i+1
    k= 0
    DO
      IF(Dim2>0) k= k +1
      j= Dims -i -k
      N= N+1
      IF(Dim2==0) EXIT
      IF(k==Dims-1 -i) EXIT
    ENDDO
    IF(i==Dims -i0) EXIT
  ENDDO
  !
  RETURN
END SUBROUTINE DiscretModel_Dim

!~ SUBROUTINE DiscretModel_Read_(vMixModel)
!~ !.called by DiscretPhase_Add in MODULE m_global_alloc
  !~ USE M_T_MixModel,ONLY: T_MixModel
  !~ !
  !~ USE M_Global_Vars,ONLY: vDiscretModel
  !~ !
  !~ TYPE(T_MixModel),INTENT(IN) :: vMixModel(:)  
  !~ !
  !~ TYPE(T_LnkDiscretModel),POINTER:: Lnk
  !~ INTEGER:: N
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_Read_begin"
  !~ !
  !~ CALL DiscretModel_BuildLnk(vMixModel,Lnk,N)
  !~ !
  !~ IF(N>0) THEN
    !~ DEALLOCATE(vDiscretModel)
    !~ ALLOCATE(vDiscretModel(N))
    !~ CALL DiscretModel_LnkToVec(Lnk,vDiscretModel)
  !~ ENDIF
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A,/)') "< DiscretModel_Read_END"
  !~ !! IF(iDebug>0) PRINT *,"DiscretModel_Read, N=", N
  !~ !
!~ ENDSUBROUTINE DiscretModel_Read_

!~ SUBROUTINE DiscretModel_BuildLnk( &
!~ !.READs the PARAMETERs for the "discretization" of a solid solution 
!~ & vMixModel, & !IN,  DATAbase of solution models
!~ & Lnk, & !POINTER
!~ & N) !OUT
  !~ USE M_IOTools
  !~ USE M_Files,     ONLY: NamFInn
  !~ USE M_T_MixModel,ONLY: T_MixModel
  !~ USE M_T_DiscretModel
  !~ !
  !~ TYPE(T_MixModel),INTENT(IN) :: vMixModel(:)
  !~ TYPE(T_LnkDiscretModel),POINTER:: Lnk
  !~ INTEGER,         INTENT(OUT):: N
  !~ !
  !~ TYPE(T_LnkDiscretModel),POINTER:: pCur
  !~ TYPE(T_DiscretModel):: DiscretModel
  !~ !
  !~ CHARACTER(LEN=512):: L,W
  !~ LOGICAL:: EoL
  !~ INTEGER:: F,ios
  !~ INTEGER:: K,I,DD
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_Read_begin"
  !~ !
  !~ N= 0
  !~ I= 0
  !~ NULLIFY(Lnk)
  !~ !
  !~ CALL GetUnit(F)
  !~ OPEN(f,FILE=TRIM(NamFInn))
  !~ DoFile: DO 
    !~ READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    !~ CALL LinToWrd(L,W,EoL)
    !~ IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    !~ CALL AppENDToEnd(L,W,EoL)
    !~ SELECT CASE(W)
    !~ !
    !~ CASE("ENDINPUT"); EXIT DoFile
    !~ !
    !~ CASE("SOLUTION.DISCRETIZE") !for READing one "block"
      !~ !! SOLUTION.DISCRETIZE
      !~ !!   NAME   BIOT
      !~ !!   NUMBER 9
      !~ !!   MODEL  BIOTITE_MG_IDEAL
      !~ !! ENDSOLUTION.DISCRETIZE
      !~ DoBlock: DO
        !~ READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        !~ CALL LinToWrd(L,W,EoL)
        !~ IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
        !~ CALL AppENDToEnd(L,W,EoL)
        !~ SELECT CASE(W)
          !~ CASE("ENDINPUT"); EXIT DoFile
          !~ CASE("END","ENDSOLUTION.DISCRETIZE")
            !~ IF(I/=0) THEN
              !~ N= N+1
              !~ !
              !~ IF(iDebug>0) WRITE(fTrc,'(A,2(A,I3))') &
              !~ & DiscretModel%Name, &
              !~ & "Dim1=  ", DiscretModel%Dim1, &
              !~ & "Model= ", DiscretModel%iMix
              !~ !
              !~ IF(N==1) THEN !_________"fill" the linked list
                !~ ALLOCATE(Lnk); NULLIFY(Lnk%next)
                !~ Lnk%Value=DiscretModel; pCur=>Lnk
              !~ ELSE
                !~ ALLOCATE(pCur%next); NULLIFY(pCur%next%next)
                !~ pCur%next%Value=DiscretModel; pCur=>pCur%next
              !~ ENDIF
              !~ !
            !~ ENDIF
            !~ EXIT DoBlock
          !~ CASE("NAME")
            !~ CALL LinToWrd(L,W,EoL); DiscretModel%Name=TRIM(W)
          !~ CASE("NUMBER")
            !~ CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,DD) 
            !~ DD= MIN(DD,MaxDiscret)
            !~ DD= MAX(DD,MinDiscret)
            !~ !
            !~ DiscretModel%Dim1=   DD -1
            !~ DiscretModel%DimTot= DiscretModel%Dim1
            !~ DiscretModel%Dim2=   0
            !~ !
            !~ IF(.not. Eol) THEN
              !~ CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,DD)
              !~ DiscretModel%Dim2= DiscretModel%Dim1
            !~ ENDIF
            !~ !
          !~ CASE("MODEL")
            !~ CALL LinToWrd(L,W,EoL)
            !~ I=0
            !~ DO K=1,SIZE(vMixModel) ! find the solution name in vMixModel%Name
              !~ IF(TRIM(W)==TRIM(vMixModel(K)%Name)) I=K
            !~ ENDDO
            !~ !
            !~ IF(vMixModel(I)%NPole<2) I=0
            !~ !
            !~ DiscretModel%iMix= I
            !~ ! IF(I==0) &
            !~ ! & CALL Stop_("SOLUTION.DISCRETIZE: "//TRIM(W)//"-> Mixing model not found ...")
        !~ ENDSELECT
      !~ ENDDO DoBlock
    !~ ENDSELECT
  !~ ENDDO DoFile
  !~ CLOSE(F)
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A,/)') "< DiscretModel_BuildLnk_END"
  !~ !
!~ ENDSUBROUTINE DiscretModel_BuildLnk

!~ SUBROUTINE DiscretModel_LnkToVec(Lnk,vDiscretModel)
  !~ USE M_T_DiscretModel,ONLY: T_DiscretModel
  !~ !
  !~ TYPE(T_LnkDiscretModel), POINTER    :: Lnk
  !~ TYPE(T_DiscretModel),    INTENT(OUT):: vDiscretModel(:)
  !~ !
  !~ TYPE(T_LnkDiscretModel),POINTER:: pCur, pPrev
  !~ INTEGER:: K
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< DiscretModel_LnkToVec_begin"
  !~ IF(iDebug>1) PRINT '(A)',"DiscretModel_LnkToVec_begin"
  !~ !
  !~ pCur=>Lnk
  !~ K=0
  !~ DO WHILE (ASSOCIATED(pCur))
    !~ K=K+1
    !~ vDiscretModel(K)=pCur%Value
    !~ pPrev=>pCur
    !~ pCur=> pCur%next
    !~ DEALLOCATE(pPrev)
    !~ !
    !~ IF(iDebug>0) WRITE(fTrc,'(I3,2A,I3)') &
    !~ & k," DiscretModel(i)%Name=",vDiscretModel(K)%Name,vDiscretModel(K)%iMix
  !~ ENDDO
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< DiscretModel_LnkToVec_END"
  !~ !
!~ ENDSUBROUTINE DiscretModel_LnkToVec

ENDMODULE M_DiscretModel_Read
