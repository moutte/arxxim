MODULE M_Test_Read
!==================================================================
! Read Block Test and Options
!==================================================================
  USE M_Trace,ONLY: iDebug,fTrc,Stop_,T_
  USE M_Kinds
  
  IMPLICIT NONE
  
  PRIVATE
  
  INTEGER :: OptComputeSize
  CHARACTER(LEN=10), ALLOCATABLE :: OptCompute(:)
  LOGICAL :: Ok_Initialized = .FALSE.
  
  PUBLIC:: Test_Read_Init
  PUBLIC:: Test_Read_OptCompute
  
CONTAINS

  FUNCTION Test_Read_OptCompute(I, Value) result (Ok)
  !================================================== 
  ! Get the COMPUTE Field number i from block TEST
  !==================================================
    IMPLICIT NONE
    LOGICAL :: Ok
    INTEGER, INTENT(IN) :: I
    CHARACTER(LEN=*), INTENT(OUT) :: Value
    !---
    IF (I>OptComputeSize) THEN
      Ok = .false.
      Value = "NONE"
    ELSE
      Ok = .true.
      Value = OptCompute(I)
    END IF
    
  END FUNCTION Test_Read_OptCompute

  !---
  
  SUBROUTINE Test_Read_Init(Ok)
  !================================================== 
  ! Init Test Block => Reinit if Necessary
  !==================================================
    IMPLICIT NONE
    LOGICAL,INTENT(OUT) :: Ok
    !--
    IF (.NOT. Ok_Initialized) CALL Test_Read_ReInit(Ok)
  
  END SUBROUTINE Test_Read_Init
  
  !---

  SUBROUTINE Test_Read_ReInit(Ok)
  !================================================== 
  ! Reinit Test block = Read Test Block
  !==================================================
    USE M_Files_Vars,ONLY:NamFInn
    USE M_IOTools 
    IMPLICIT NONE
    !-----
    LOGICAL,INTENT(OUT) :: Ok
    !-----
    LOGICAL :: sEOL,EoL
    INTEGER :: nCompute
    CHARACTER(LEN=512):: L,W
    INTEGER :: ios,F
    INTEGER :: i
    !-----
    Ok = .FALSE.
    nCompute =0

    IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< Test_ReadInput"

    !// Estimate the maximal number of values to store
    CALL GetUnit(f)
    OPEN(f,FILE=TRIM(NamFInn),STATUS='OLD')

    nCompute = 0
    DO 
      READ(f,'(A)',IOSTAT=ios) L
      IF(ios/=0) EXIT
      CALL LinToWrd(L,W,EoL)
      IF(W(1:1)=='!') CYCLE
      IF(TRIM(W)=="ENDINPUT") EXIT
      IF(TRIM(W)=="COMPUTE") nCompute = nCompute +1
    ENDDO
    CLOSE(f)
    
    !~ print *,'nCompute=',nCompute
    !~ print *,'..DEBBUGG'  ;  pause
    
    !// Allocate the table
    IF(ALLOCATED(OptCompute)) DEALLOCATE(OptCompute)
    ALLOCATE(OptCompute(10))

    !// Read the File and store the values
    CALL GetUnit(f)
    OPEN(f,FILE=TRIM(NamFInn),STATUS='OLD')
    !
    nCompute = 0
    !
    DoFile: DO
      !
      READ(f,'(A)',IOSTAT=ios) L
      IF(ios/=0) EXIT DoFile
      CALL LinToWrd(L,W,EoL)
      IF(W(1:1)=='!') CYCLE DoFile
      CALL AppendToEnd(L,W,EoL)
      
      SELECT CASE(TRIM(W))
          
      CASE("ENDINPUT") ; EXIT DoFile
      
      CASE("TEST")
        Ok= .TRUE.
        !
        DoTest: DO
          !
          READ(F,'(A)',IOSTAT=ios) L
          IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoTest
          CALL AppendToEnd(L,W,EoL)
          
          SELECT CASE(TRIM(W))
          
          CASE("ENDINPUT")
            EXIT DoFile
          CASE("END","ENDTEST")
            EXIT DoTest
          
          CASE("COMPUTE")
            CALL LinToWrd(L,W,sEol) 
            nCompute = nCompute + 1
            OptCompute(nCompute)= TRIM(W)
          
          CASE DEFAULT
            Ok= .FALSE.
            PRINT *,TRIM(W)//" =INVALID keyword in TEST block"
          
          END SELECT
          !
        ENDDO Dotest
        !
      END SELECT
      !
    ENDDO DoFile
    CLOSE(f)
    !
    OptComputeSize = nCompute
    Ok_Initialized = .true.
    
    do nCompute=1,OptComputeSize
      WRITE(fTrc,'(2X,A)') OptCompute(nCompute)
    enddo
    !
    IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "</ Test_ReadInput"

  ENDSUBROUTINE Test_Read_ReInit

END MODULE M_Test_Read
