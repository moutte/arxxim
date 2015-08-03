MODULE M_TPCond_Read
  USE M_Kinds
  USE M_Trace

  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC:: TPpath_Read
  PUBLIC:: TPgrid_Read
  PUBLIC:: TPgrid_Build
  !
CONTAINS

SUBROUTINE TPpath_Read(TdgK0,Pbar0)
!--
!-- scan the input file, reads block TP.TABLE, if present
!-- -> allocate & build vTPpath
!--
  USE M_IoTools !,ONLY:GetUnit,dimV
  USE M_Files,     ONLY: NamFInn
  !
  USE M_Dtb_Const, ONLY: T_CK,Tref,Pref
  USE M_Dtb_Calc,  ONLY: Dtb_TP_Check
  USE M_T_Tpcond,  ONLY: T_TPCond
  USE M_Fluid_Calc,ONLY: Eos_H2O_Psat
  !
  USE M_Dtb_Vars,  ONLY: DtbFormat,DtbLogK_vTPCond,Psat_Auto
  USE M_Path_Vars, ONLY: vTPpath,DimPath
  !
  REAL(dp),INTENT(IN),OPTIONAL:: TdgK0,Pbar0
  != default values when no TP.TABLE or TP.PATH is found
  !
  LOGICAL :: EoL,Ok
  INTEGER :: i,mDum,ios,f,N !,m1,m2
  REAL(dp):: TdgK,Pbar
  CHARACTER(LEN=512):: L,W,W1
  CHARACTER(LEN=80) :: Msg
  !
  REAL(dp):: vX(dimV)
  TYPE(T_TPCond):: vCond(dimV)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< TPpath_Read"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  Ok=.FALSE.
  !
  vCond(:)%Name= "x"
  vCond(:)%Pbar= Pref
  !
  DoFile: DO
    !
    READ(f,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT  DoFile
    !
    CASE("TP.TABLE","TPTABLE","TP.PATH")

      IF(iDebug>0 .and. &
      &   DtbFormat=="LOGKTBL" .and. &
      &  (TRIM(W)=="TP.TABLE" .or. TRIM(W)=="TPTABLE")) THEN
        CALL Warning_("with logK table, USE TP.PATH, to distinguish from TP.TABLE")
      ENDIF

      Ok=  .TRUE.
      N= dimV
      DoTPtable: DO

        READ(f,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,Eol)
        IF(W(1:1)=='!') CYCLE DoTPtable
        !WRITE(fTrc,'(A)') L
        CALL AppendToEnd(L,W,EoL)

        SELECT CASE(W)

          !! CASE("ENDINPUT"); EXIT DoFile
          CASE("END","ENDTPTABLE","ENDTP.TABLE","ENDTP.PATH")
            EXIT DoTPtable
          !---------------!

          CASE("TDGC", "TDGK", "PBAR","PMPA", "TDEGC","TDEGK")
            !
            !! CALL ReadRValsV(L,mDum,vX)
            !
            i=0
            DO
              CALL LinToWrd(L,W1,EoL)
              i=i+1
              IF(i>DimV) EXIT
              IF(TRIM(W1)=="PSAT") THEN
                IF(vCond(i)%TdgC>Zero) THEN
                  CALL Eos_H2O_Psat(vCond(i)%TdgC+T_CK,vX(i))
                  !PRINT *,"TdgC= ", vCond(i)%TdgC
                  !PRINT *,"psat= ", vX(i)  ;  CALL PaUSE_
                ELSE
                  CALL Stop_("error in TP.TABLE: invalid TdgC for Psat computation")
                ENDIF
              ELSE
                CALL WrdToReal(W1,vX(i))
              ENDIF
              IF(EoL) EXIT
            ENDDO
            !
            mDum= i
            N=MIN(N,mDum)
            IF(iDebug>0) &
            & WRITE(fTrc,'(A,A1,A,2I3)') TRIM(W),T_," DIM= ", mDum, N
            !
            SELECT CASE(TRIM(W))
            !
            CASE("TDGC");  vCond(:)%TdgC= vX(:)
            CASE("TDGK");  vCond(:)%TdgC= vX(:)-T_CK
            !
            CASE("PBAR");  vCond(:)%Pbar= vX(:)       !pressure in bar
            CASE("PMPA");  vCond(:)%Pbar= vX(:)*1.D-1 !pressure in MegaPascal
            !
            CASE("TDEGC"); vCond(:)%TdgC=vX(:)
              CALL Warning_("Depreciated Keyword TDEGC")
            CASE("TDEGK"); vCond(:)%TdgC=vX(:)-T_CK
              CALL Warning_("Depreciated Keyword TDEGK")
            !
            END SELECT
          !---------------!

          CASE("NAME")
            i=0
            DoName: DO
              CALL LinToWrd(L,W,EoL)
              i=i+1
              IF(I>DimV) EXIT DoName
              vCond(i)%Name=TRIM(W)
              IF(EoL) EXIT DoName
            ENDDO DoName
          !---------------!

          CASE DEFAULT
            CALL Stop_(TRIM(W)//"<<unknown Keyword...")

        END SELECT !ENDIF
        !ENDIF
      ENDDO DoTPtable
    END SELECT
    !
  ENDDO DoFile
  CLOSE(f)
  !
  IF(Ok) THEN

    DO I=1,N
      TdgK= vCond(I)%TdgC +T_CK
      Pbar= vCond(I)%Pbar
      !
      CALL Dtb_TP_Check(DtbFormat,DtbLogK_vTPCond,Psat_Auto,TdgK,Pbar,Ok,Msg)
      !
      IF(.NOT. Ok) CALL Stop_(TRIM(Msg))
      !
      vCond(I)%TdgC= TdgK -T_CK
      vCond(I)%Pbar= Pbar
      !
      IF(iDebug>0) &
      & WRITE(fTrc,'(2(G12.3,1X))') vCond(I)%TdgC, vCond(I)%Pbar !!i,vCond(i)%Name
    ENDDO

  ELSE

    N=1
    IF(PRESENT(TdgK0)) THEN ;  vCond(1)%TdgC= TdgK0 -T_CK
    ELSE                    ;  vCond(1)%TdgC= Tref  -T_CK
    ENDIF
    IF(PRESENT(Pbar0)) THEN ;  vCond(1)%Pbar= Pbar0
    ELSE                    ;  vCond(1)%Pbar= Pref
    ENDIF

  ENDIF
  !
  IF(ALLOCATED(vTPpath)) DEALLOCATE(vTPpath)
  ALLOCATE(vTPpath(1:N))
  vTPpath(1:N)= vCond(1:N)
  !
  DimPath= N
  !
  !!IF(.NOT.Ok) CALL Stop_("Block TP.PATH not Found ...!!!")
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ TPpath_Read"
  !
ENDSUBROUTINE TPpath_Read

SUBROUTINE TPgrid_Read( &
& Ok, &
& T_Min,T_Max,T_ratio,T_delta, &
& P_Min,P_Max,P_ratio,P_delta)
  !
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  USE M_Dtb_Const,ONLY: T_CK
  USE M_Files,    ONLY: NamFInn
  !
  LOGICAL, INTENT(OUT):: Ok
  REAL(dp),INTENT(OUT):: T_Min,T_Max,T_ratio,T_delta,P_Min,P_Max,P_ratio,P_delta
  !
  CHARACTER(LEN=512):: L,W,W1
  LOGICAL :: EoL
  INTEGER :: F,ios
  REAL(dp):: rBegin,rFinal,rRatio,rDelta
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dtb_Read_TPgrid"
  !
  Ok= .FALSE.
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))
  !
  DoFile: DO

    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines

    CALL AppENDToEnd(L,W,EoL)
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT DoFile
    !
    CASE("TP.GRID","TPGRID") !!!!!!canevas for reading one "block"
      !... I=0
      Ok= .TRUE.
      !
      DoBlock: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT"); EXIT DoFile
        !
        CASE("END","ENDTP.GRID","ENDTPGRID"); EXIT DoBlock
        !
        CASE DEFAULT
          CALL Pause_("< WARNING!!! "//TRIM(W)//"=UNKNOWN keyword in TPgrid")
        !
        CASE("TDGC","TDGK","PBAR") !,"PMPA")
          !!PRINT '(A)',TRIM(L)  ; PAUSE_
          rBegin=Zero; rFinal=Zero; rRatio=One; rDelta=Zero
          !
          DoLine: DO
            IF(EoL) EXIT DoLine
            CALL LinToWrd(L,W1,EoL)
            SELECT CASE(W1)
              CASE("INITIAL"); CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rBegin)
              CASE("FINAL");   CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rFinal)
              CASE("RATIO");   CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rRatio)
              CASE("DELTA");   CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rDelta)
              CASE DEFAULT;    CALL Stop_(W1//"= unknown keyword in TPgrid !!!") !STOP
            END SELECT
          ENDDO DoLine
          !
          SELECT CASE(W)
          CASE("TDGC")
            T_Min=   rBegin ;      T_Max=   rFinal
            T_ratio= rRatio ;      T_delta= rDelta
          CASE("TDGK")
            T_Min=   rBegin-T_CK;  T_Max=   rFinal-T_CK
            T_ratio= rRatio ;      T_delta= rDelta
          CASE("PBAR")
            P_Min=   rBegin ;      P_Max=   rFinal
            P_ratio= rRatio ;      P_delta= rDelta
          END SELECT
        !ENDCASE("TDGC","PBAR","PMPA")
        !
        END SELECT
        !
      ENDDO DoBlock
    !
    END SELECT
  ENDDO DoFile
  CLOSE(F)

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dtb_Read_TPgrid"

  RETURN
ENDSUBROUTINE TPgrid_Read

SUBROUTINE TPgrid_Build(Ok)
  USE M_T_TPcond,   ONLY: T_TPCond
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_Fluid_Calc, ONLY: Eos_H2O_Psat
  !
  USE M_Path_Vars,  ONLY: vTPpath
  !
  LOGICAL,INTENT(OUT):: Ok
  !
  REAL(dp):: T_Min, T_Max, T_delta, T_ratio
  REAL(dp):: P_Min, P_Max, P_delta, P_ratio
  REAL(dp):: TdgC,TdgK,Pbar,Psat
  !
  INTEGER:: N
  !
  CALL TPgrid_Read( Ok, &
  & T_Min,T_Max,T_ratio,T_delta, &
  & P_Min,P_Max,P_ratio,P_delta)
  !
  ! when a specific TP.GRID block is not found,
  ! use a simple TP.TABLE
  IF(.NOT. Ok) THEN
    != IF TP.GRID block not found
    !N= SIZE(vTPcond)
    !IF(ALLOCATED(vTPgrid)) DEALLOCATE(vTPgrid)
    !ALLOCATE(vTPgrid(1:N))
    !vTPgrid= vTPcond
    RETURN
  ENDIF
  !
  IF(iDebug==4) THEN
    PRINT *,"TPgrid_Build"
    PRINT '(A,4G15.6)',"T_Min,T_Max,T_ratio,T_delta",T_Min,T_Max,T_ratio,T_delta
    PRINT '(A,4G15.6)',"P_Min,P_Max,P_ratio,P_delta",P_Min,P_Max,P_ratio,P_delta
    CALL Pause_
  ENDIF

  !-------------------------------- first, count number of T,P points --
  N= 0
  TdgC= T_Min
  DO
    Pbar= P_Min
    DO
      TdgK= TdgC +T_CK
      Psat= Zero
      IF(TdgK<=647.25D0) CALL Eos_H2O_Psat(TdgK,Psat)
      IF(Pbar>Psat) N= N+1
      IF(P_ratio>One) THEN ; Pbar= Pbar * P_ratio
      ELSE                 ; Pbar= Pbar + P_delta
      ENDIF
      IF(Pbar>P_Max) EXIT
    ENDDO
    TdgC= TdgC + T_delta
    IF(TdgC>T_max) EXIT
  ENDDO
  !
  IF(ALLOCATED(vTPpath)) DEALLOCATE(vTPpath)
  ALLOCATE(vTPpath(1:N))

  !----------------------------------------------- then, fill vTPpath --
  TdgC= T_Min
  N= 0
  DO
    Pbar= P_Min
    DO
      TdgK= TdgC +T_CK
      !PRINT '(F7.2,1X,F7.2)',TdgC,Pbar
      Psat= Zero
      IF (TdgK<=647.25D0) CALL Eos_H2O_Psat(TdgK,Psat)
      IF(Pbar>Psat) THEN
        N= N+1
        vTPpath(N)%TdgC= TdgC
        vTPpath(N)%Pbar= Pbar
      ENDIF
      IF(P_ratio>One) THEN ; Pbar= Pbar * P_ratio
      ELSE                 ; Pbar= Pbar + P_delta
      ENDIF
      IF(Pbar>P_Max) EXIT
    ENDDO
    TdgC= TdgC + T_delta
    IF(TdgC>T_max) EXIT
  ENDDO

  RETURN
ENDSUBROUTINE TPgrid_Build

ENDMODULE M_TPcond_Read
