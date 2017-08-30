module M_GETformat
  !//-------------------------------------------------------------------
  !// Utils to obtain the format associated to a string of a number
  !// Exemple.
  !// call getformat ('1.23456D+78',  myformat)
  !// Resultat : myformat = (1PD11.5)
  !//-------------------------------------------------------------------

  implicit none
  private

  integer, parameter :: KINT = 1 ! integer
  integer, parameter :: KFIX = 2 ! fixed point real
  integer, parameter :: KEXP = 3 ! exponent type real
  integer, parameter :: KDBL = 4 ! exponent type double


  private :: isnum
  private :: obtfmt
  private :: NBRCHF

  interface getformat
    module procedure obtfmt
  end interface

  public :: getformat
  public :: tstfmt

contains

  subroutine tstfmt
    implicit none
    character (len=64) :: ztxt
    character (len=64) :: zfmt

    write(*,*) "INPUT STRING"
    read (*, "(A)") ztxt
    call obtfmt (ztxt, zfmt)
    write (*, *) zfmt

    stop "PERFECT"
  end subroutine tstfmt

  integer function isnum (ZVAL)
    !  Verify that a character string represents a numerical value
    implicit none
    integer :: NUM, NMTS, NEXP, KMTS, ifEXP, ICHR
    character (Len=*), intent (in) :: ZVAL
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
    ifEXP = 0
    !
    ! loop over characters
    !
    ICHR = 0
    do
      if (ICHR >= len(ZVAL)) then
        !
        ! last check
        !
        if (NMTS == 0) exit
        if (NUM >= KEXP .and. NEXP == 0) exit
        isnum = NUM
        return
      end if
      ICHR = ICHR + 1
      select case (ZVAL(ICHR:ICHR))
        !
        ! process blanks
        !
        case (' ')
          continue
          !
          ! process digits
          !
        case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
          if (NUM == 0) NUM = KINT
          if (NUM < KEXP) then
            NMTS = NMTS + 1
          else
            NEXP = NEXP + 1
          end if
          !
          ! process signs
          !
        case ('+', '-')
          if (NUM == 0) then
            if (KMTS > 0) exit
            KMTS = 1
            NUM = KINT
          else
            if (NUM < KEXP) exit
            if (ifEXP > 0) exit
            ifEXP = 1
          end if
          !
          ! process decimal point
          !
        case ('.')
          if (NUM /= KINT .and. ICHR /= 1) exit
          NUM = KFIX
          !
          ! process exponent
          !
        case ('e', 'E')
          if (NUM >= KEXP) exit
          if (NMTS == 0) exit
          NUM = KEXP
          !
        case ('d', 'D')
          if (NUM >= KEXP) exit
          if (NMTS == 0) exit
          NUM = KDBL
          !
          ! any other character means the string is non-numeric
          !
        case default
          exit
          !
      end select
    end do
    !
    ! if this point is reached, the string is non-numeric
    !
    isnum = 0
    return
  end function isnum

  subroutine obtfmt (ZVAL, zfmt)
    implicit none
    integer :: LVAL, LFMT, KASE, NCHR , NCHF,  NCHR1,  NCHFF,  NCHR2,  IPNT, NCHFP,  NCHR0
    !  Find out what Fortran format was used to write a numerical value
    character (Len=*), intent (in) :: ZVAL ! the string
    character (Len=*), intent (out) :: zfmt ! the format
    ! ____________________________________________________
    character (Len=1) :: ZNUM0, ZNUM1, ZNUM2 ! to write the numbers
    ! of digits of the integers
    ! used in the format
    character (Len=27) :: ZFMTW ! The format to write the format ...

    !
    ! initialise
    !
    LVAL = len (ZVAL)
    LFMT = len (zfmt)
    !
    !  Big switching place
    !
    KASE = isnum (ZVAL)
    select case (KASE)
      !
      ! non numeric (A Format)
      ! ____________________________________________________
      !
      case (0)
        NCHR = LVAL
        NCHR1 = NBRCHF (NCHR)
        !
        !    The format is (Axxx), we will write it with ('(A',Iw,')')
        !    Lets create the latter format, ZFMTW, with w = NCHR1
        !
        if (NCHR1+3 > LFMT) then
          write (*, *) "Argument string ZFMT too short"
        else
          if (NCHR1 > 0 .and. NCHR1 < 10) then
            write (ZNUM1, "(I1)") NCHR1
            ZFMTW = "('(A',I" // ZNUM1 // ",')')"
            write (zfmt, ZFMTW) NCHR
          else
            write (*, *) "Doesn't a string length of",&
              & NCHR, " seem strange ?"
          end if
        end if
        !
        ! integer
        ! ____________________________________________________
        !
      case (KINT)
        NCHF = len_trim (ZVAL)
        !
        ! If it looks too long, remove leading blanks
        !
        if (NCHF > 20) then
          NCHF = len_trim (adjustl(ZVAL))
        end if
        !
        NCHR1 = NBRCHF (NCHF)
        !
        !    The format is (Ixxx), we will write it with  ('(I',Iw,')')
        !    Lets create the latter format, ZFMTW, with w = NCHR1
        !
        if (NCHR1+3 > LFMT) then
          write (*, *) "Argument string ZFMT too short"
        else
          if (NCHR1 > 0 .and. NCHR1 < 10) then
            write (ZNUM1, "(I1)") NCHR1
            ZFMTW = "('(I',I" // ZNUM1 // ",')')"
            write (zfmt, ZFMTW) NCHF
          else
            write (*, *) "isn't an integer of ", NCHF, &
              & "digits a strange idea ?"
          end if
        end if
        !
        ! real, fixed point form
        ! ____________________________________________________
        !
      case (KFIX)
        NCHF = len_trim (ZVAL)
        NCHFF = NCHF - INDEX (ZVAL, '.')
        !
        ! If it looks too long, remove leading blanks
        !
        if (NCHF > 25) then
          NCHF = len_trim (adjustl(ZVAL))
        end if
        !
        NCHR1 = NBRCHF (NCHF)
        NCHR2 = NBRCHF (NCHFF)
        !
        !    The format is (Fxx.yy), we will write it with  ('(F',Iw,'.',Id,')')
        !    Lets create the latter format, ZFMTW, with w = NCHR1
        !    and d = NCHR2, obtained from the position of the decimal point
        !
        if (NCHR1+NCHR2+4 > LFMT) then
          write (*, *) "Argument string ZFMT too short"
        else
          if (NCHR1 > 0 .and. NCHR1 < 10) then
            write (ZNUM1, "(I1)") NCHR1
            write (ZNUM2, "(I1)") NCHR2
            ZFMTW = "('(F',I" // ZNUM1 // ",'.',I" // ZNUM2 //&
              & ",')')"
            write (zfmt, ZFMTW) NCHF, NCHFF
          else
            write (*, *) "isn't a real of ", NCHF, &
            & "digits a strange idea ?"
          end if
        end if
        !
        ! real, exponent form
        ! ____________________________________________________
        !
      case (KEXP, KDBL)
        NCHF = len_trim (ZVAL)
        if (KASE == 3) then
          NCHFF = Max (INDEX(ZVAL, 'E'), INDEX(ZVAL, 'e')) - 1 -&
            & INDEX (ZVAL, '.')
        else
          NCHFF = Max (INDEX(ZVAL, 'D'), INDEX(ZVAL, 'd')) - 1 -&
            & INDEX (ZVAL, '.')
        end if
        IPNT = INDEX (ZVAL, '.')
        if (IPNT > 0) then
          NCHFP = IPNT - VERifY (ZVAL, " +-")
        else
          NCHFP = NCHFF
          NCHFF = 0
        end if
        !
        ! If it looks too long, remove leading blanks
        !
        if (NCHF > 30) then
          NCHF = len_trim (adjustl(ZVAL))
        end if
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
        if (NCHF+5 > LFMT) then
          write (*, *) "Argument string ZFMT too short"
        else
          if (NCHR1 > 0 .and. NCHR1 < 10) then
            write (ZNUM0, "(I1)") NCHR0
            write (ZNUM1, "(I1)") NCHR1
            write (ZNUM2, "(I1)") NCHR2
            if (KASE == 3) then
              ZFMTW = "('(',I" // ZNUM0 // ",'PE',I" // ZNUM1 //&
                & ",'.',I" // ZNUM2 // ",')')"
            else
              ZFMTW = "('(',I" // ZNUM0 // ",'PD',I" // ZNUM1 //&
                & ",'.',I" // ZNUM2 // ",')')"
            end if
            write (zfmt, ZFMTW) NCHFP, NCHF, NCHFF
          else
            write (*, *) "isn't a real of ", NCHF, &
              & "digits a strange idea ?"
          end if
        end if
        !
        !
      case default
        write (*, *) "Type ", KASE, " not known"
        !
    end select
    return
  end subroutine obtfmt

  integer function NBRCHF (JVAL)
    !  Number of characters (digits and minus sign) to display JVAL
    implicit none
    integer ::  NCHF,  JVALA
    integer, intent (in) :: JVAL ! the value
    ! ____________________________________________________
    !
    ! Compute with integers to avoid precision problems
    ! with logarithms
    !
    ! Start with 1, [+1 when negative]
    ! ____________________________________________________
    !
    if (JVAL < 0) then
      NCHF = 2
      JVALA = - JVAL
    else
      NCHF = 1
      JVALA = JVAL
    end if
    !
    !         + 1 per overpassing of power of 10
    !
    do
      if (JVALA < 10) exit
      JVALA = JVALA / 10
      NCHF = NCHF + 1
    end do
    NBRCHF = NCHF
    return
  end function NBRCHF

end module M_GETformat
