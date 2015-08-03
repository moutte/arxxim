MODULE M_GETFORMAT
  !//-------------------------------------------------------------------
  !// Utils to obtain the format associated to a string of a number
  !// Exemple.
  !// CALL getformat ('1.23456D+78',  myformat)
  !// Resultat : myformat = (1PD11.5)
  !//-------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: KINT = 1 ! integer
  INTEGER, PARAMETER :: KFIX = 2 ! fixed point real
  INTEGER, PARAMETER :: KEXP = 3 ! exponent type real
  INTEGER, PARAMETER :: KDBL = 4 ! exponent type double


  PRIVATE :: isnum
  PRIVATE :: obtfmt
  PRIVATE :: NBRCHF

  INTERFACE getformat
    module procedure obtfmt
  END INTERFACE

  PUBLIC :: getformat
  PUBLIC :: tstfmt

CONTAINS

  SUBROUTINE tstfmt
    IMPLICIT NONE
    CHARACTER (LEN=64) :: ztxt
    CHARACTER (LEN=64) :: zfmt

    WRITE(*,*) "INPUT STRING"
    READ (*, "(A)") ztxt
    CALL obtfmt (ztxt, zfmt)
    WRITE (*, *) zfmt

    STOP "PERFECT"
  END SUBROUTINE tstfmt

  INTEGER FUNCTION isnum (ZVAL)
    !  Verify that a character string represents a numerical value
    IMPLICIT NONE
    INTEGER :: NUM, NMTS, NEXP, KMTS, IFEXP, ICHR
    CHARACTER (Len=*), INTENT (IN) :: ZVAL
    ! __________________________________________________
    !  return : 0 = non-numeric string
    !        else = code as defined in module codnum
    ! __________________________________________________
    !
    ! initialise
    !
    NUM = 0
    NMTS = 0
    NEXP = 0
    KMTS = 0
    IFEXP = 0
    !
    ! loop over characters
    !
    ICHR = 0
    DO
      IF (ICHR >= LEN(ZVAL)) THEN
        !
        ! last check
        !
        IF (NMTS == 0) EXIT
        IF (NUM >= KEXP .AND. NEXP == 0) EXIT
        isnum = NUM
        RETURN
      END IF
      ICHR = ICHR + 1
      SELECT CASE (ZVAL(ICHR:ICHR))
        !
        ! process blanks
        !
        CASE (' ')
          CONTINUE
          !
          ! process digits
          !
        CASE ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
          IF (NUM == 0) NUM = KINT
          IF (NUM < KEXP) THEN
            NMTS = NMTS + 1
          ELSE
            NEXP = NEXP + 1
          END IF
          !
          ! process signs
          !
        CASE ('+', '-')
          IF (NUM == 0) THEN
            IF (KMTS > 0) EXIT
            KMTS = 1
            NUM = KINT
          ELSE
            IF (NUM < KEXP) EXIT
            IF (IFEXP > 0) EXIT
            IFEXP = 1
          END IF
          !
          ! process decimal point
          !
        CASE ('.')
          IF (NUM /= KINT .AND. ICHR /= 1) EXIT
          NUM = KFIX
          !
          ! process exponent
          !
        CASE ('e', 'E')
          IF (NUM >= KEXP) EXIT
          IF (NMTS == 0) EXIT
          NUM = KEXP
          !
        CASE ('d', 'D')
          IF (NUM >= KEXP) EXIT
          IF (NMTS == 0) EXIT
          NUM = KDBL
          !
          ! any other character means the string is non-numeric
          !
        CASE DEFAULT
          EXIT
          !
      END SELECT
    END DO
    !
    ! if this point is reached, the string is non-numeric
    !
    isnum = 0
    RETURN
  END FUNCTION isnum

  SUBROUTINE obtfmt (ZVAL, zfmt)
    IMPLICIT NONE
    INTEGER :: LVAL, LFMT, KASE, NCHR , NCHF,  NCHR1,  NCHFF,  NCHR2,  IPNT, NCHFP,  NCHR0
    !  Find out what Fortran format was used to write a numerical value
    CHARACTER (Len=*), INTENT (IN) :: ZVAL ! the string
    CHARACTER (Len=*), INTENT (OUT) :: zfmt ! the format
    ! ____________________________________________________
    CHARACTER (Len=1) :: ZNUM0, ZNUM1, ZNUM2 ! to write the numbers
    ! of digits of the integers
    ! used in the format
    CHARACTER (Len=27) :: ZFMTW ! The format to write the format ...

    !
    ! initialise
    !
    LVAL = LEN (ZVAL)
    LFMT = LEN (zfmt)
    !
    !  Big switching place
    !
    KASE = isnum (ZVAL)
    SELECT CASE (KASE)
      !
      ! non numeric (A Format)
      ! ____________________________________________________
      !
      CASE (0)
        NCHR = LVAL
        NCHR1 = NBRCHF (NCHR)
        !
        !    The format is (Axxx), we will write it with ('(A',Iw,')')
        !    Lets create the latter format, ZFMTW, with w = NCHR1
        !
        IF (NCHR1+3 > LFMT) THEN
          WRITE (*, *) "Argument string ZFMT too short"
        ELSE
          IF (NCHR1 > 0 .AND. NCHR1 < 10) THEN
            WRITE (ZNUM1, "(I1)") NCHR1
            ZFMTW = "('(A',I" // ZNUM1 // ",')')"
            WRITE (zfmt, ZFMTW) NCHR
          ELSE
            WRITE (*, *) "Doesn't a string length of",&
              & NCHR, " seem strange ?"
          END IF
        END IF
        !
        ! integer
        ! ____________________________________________________
        !
      CASE (KINT)
        NCHF = LEN_TRIM (ZVAL)
        !
        ! If it looks too long, remove leading blanks
        !
        IF (NCHF > 20) THEN
          NCHF = LEN_TRIM (ADJUSTL(ZVAL))
        END IF
        !
        NCHR1 = NBRCHF (NCHF)
        !
        !    The format is (Ixxx), we will write it with  ('(I',Iw,')')
        !    Lets create the latter format, ZFMTW, with w = NCHR1
        !
        IF (NCHR1+3 > LFMT) THEN
          WRITE (*, *) "Argument string ZFMT too short"
        ELSE
          IF (NCHR1 > 0 .AND. NCHR1 < 10) THEN
            WRITE (ZNUM1, "(I1)") NCHR1
            ZFMTW = "('(I',I" // ZNUM1 // ",')')"
            WRITE (zfmt, ZFMTW) NCHF
          ELSE
            WRITE (*, *) "isn't an integer of ", NCHF, &
              & "digits a strange idea ?"
          END IF
        END IF
        !
        ! real, fixed point form
        ! ____________________________________________________
        !
      CASE (KFIX)
        NCHF = LEN_TRIM (ZVAL)
        NCHFF = NCHF - INDEX (ZVAL, '.')
        !
        ! If it looks too long, remove leading blanks
        !
        IF (NCHF > 25) THEN
          NCHF = LEN_TRIM (ADJUSTL(ZVAL))
        END IF
        !
        NCHR1 = NBRCHF (NCHF)
        NCHR2 = NBRCHF (NCHFF)
        !
        !    The format is (Fxx.yy), we will write it with  ('(F',Iw,'.',Id,')')
        !    Lets create the latter format, ZFMTW, with w = NCHR1
        !    and d = NCHR2, obtained from the position of the decimal point
        !
        IF (NCHR1+NCHR2+4 > LFMT) THEN
          WRITE (*, *) "Argument string ZFMT too short"
        ELSE
          IF (NCHR1 > 0 .AND. NCHR1 < 10) THEN
            WRITE (ZNUM1, "(I1)") NCHR1
            WRITE (ZNUM2, "(I1)") NCHR2
            ZFMTW = "('(F',I" // ZNUM1 // ",'.',I" // ZNUM2 //&
              & ",')')"
            WRITE (zfmt, ZFMTW) NCHF, NCHFF
          ELSE
            WRITE (*, *) "isn't a real of ", NCHF, &
            & "digits a strange idea ?"
          END IF
        END IF
        !
        ! real, exponent form
        ! ____________________________________________________
        !
      CASE (KEXP, KDBL)
        NCHF = LEN_TRIM (ZVAL)
        IF (KASE == 3) THEN
          NCHFF = Max (INDEX(ZVAL, 'E'), INDEX(ZVAL, 'e')) - 1 -&
            & INDEX (ZVAL, '.')
        ELSE
          NCHFF = Max (INDEX(ZVAL, 'D'), INDEX(ZVAL, 'd')) - 1 -&
            & INDEX (ZVAL, '.')
        END IF
        IPNT = INDEX (ZVAL, '.')
        IF (IPNT > 0) THEN
          NCHFP = IPNT - VERIFY (ZVAL, " +-")
        ELSE
          NCHFP = NCHFF
          NCHFF = 0
        END IF
        !
        ! If it looks too long, remove leading blanks
        !
        IF (NCHF > 30) THEN
          NCHF = LEN_TRIM (ADJUSTL(ZVAL))
        END IF
        !
        NCHR0 = NBRCHF (NCHFP)
        NCHR1 = NBRCHF (NCHF)
        NCHR2 = NBRCHF (NCHFF)
        !
        !    The chosen format is (zPExx.yy), we will write it with
        !     ('(',Ik,'PE',Iw,'.',Id,')')
        !    Lets create the latter format, ZFMTW, with
        !    k= NCHR0, w = NCHR1, d = NCHR2
        !
        IF (NCHF+5 > LFMT) THEN
          WRITE (*, *) "Argument string ZFMT too short"
        ELSE
          IF (NCHR1 > 0 .AND. NCHR1 < 10) THEN
            WRITE (ZNUM0, "(I1)") NCHR0
            WRITE (ZNUM1, "(I1)") NCHR1
            WRITE (ZNUM2, "(I1)") NCHR2
            IF (KASE == 3) THEN
              ZFMTW = "('(',I" // ZNUM0 // ",'PE',I" // ZNUM1 //&
                & ",'.',I" // ZNUM2 // ",')')"
            ELSE
              ZFMTW = "('(',I" // ZNUM0 // ",'PD',I" // ZNUM1 //&
                & ",'.',I" // ZNUM2 // ",')')"
            END IF
            WRITE (zfmt, ZFMTW) NCHFP, NCHF, NCHFF
          ELSE
            WRITE (*, *) "isn't a real of ", NCHF, &
              & "digits a strange idea ?"
          END IF
        END IF
        !
        !
      CASE DEFAULT
        WRITE (*, *) "Type ", KASE, " not known"
        !
    END SELECT
    RETURN
  END SUBROUTINE obtfmt

  INTEGER FUNCTION NBRCHF (JVAL)
    !  Number of characters (digits and minus sign) to display JVAL
    IMPLICIT NONE
    INTEGER ::  NCHF,  JVALA
    INTEGER, INTENT (IN) :: JVAL ! the value
    ! ____________________________________________________
    !
    ! Compute with integers to avoid precision problems
    ! with logarithms
    !
    ! Start with 1, [+1 when negative]
    ! ____________________________________________________
    !
    IF (JVAL < 0) THEN
      NCHF = 2
      JVALA = - JVAL
    ELSE
      NCHF = 1
      JVALA = JVAL
    END IF
    !
    !         + 1 per overpassing of power of 10
    !
    DO
      IF (JVALA < 10) EXIT
      JVALA = JVALA / 10
      NCHF = NCHF + 1
    END DO
    NBRCHF = NCHF
    RETURN
  END FUNCTION NBRCHF

END MODULE M_GETFORMAT
