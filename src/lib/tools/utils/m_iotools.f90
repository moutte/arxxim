!.module for basic output procedures
MODULE M_IOTools
  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: cUpper, cLower
  PUBLIC:: Str_Append
  PUBLIC:: Str_Upper
  PUBLIC:: CarToInt,WrdToInt
  PUBLIC:: WrdToReal,ReadRValsV
  PUBLIC:: IntToStr3,IntToStr4,FIntToStr3,FIntToStr4
  PUBLIC:: In_Int,In_Car,In_Str
  PUBLIC:: LinToWrd
  PUBLIC:: AppendToEnd
  PUBLIC:: OutStrVec
  PUBLIC:: GetUnit
  PUBLIC:: OpenFile
  PUBLIC:: CloseFile
  !
  INTEGER,  PARAMETER,PUBLIC:: dimV=32 !max.dimension for routine ReadRValsV
  !
  CHARACTER,PARAMETER:: T_= CHAR(9)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ModIORead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ModIORead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ModIORead

character function cUpper(c)
  character,intent(in):: c
  integer:: j,d
  d= ichar('A') - ichar('a')
  cUpper= c
  j= ichar(c)
  if(j>=ichar('a') .AND. j<=ichar('z')) cUpper= char(j+d)
end function cUpper

character function cLower(c)
  character,intent(in):: c
  integer:: j,d
  d= ichar('A') - ichar('a')
  cLower= c
  j= ichar(c)
  if(j>=ichar('A') .AND. j<=ichar('Z')) cLower= char(j-d)
end function cLower

SUBROUTINE GetUnit(F) !returns a free unit number.
!Author:John Burkardt
!A "free" FORTRAN unit number is an integer between 1 and 99
!which is not currently associated with an I/O device.
!A free FORTRAN unit number is needed in order to open a file with the OPEN command.
!IUNIT=0
!no free FORTRAN unit could be found, although all 99 units were checked
  IMPLICIT NONE
  INTEGER,INTENT(OUT)::F
  INTEGER:: i, ios
  LOGICAL:: lOpen
  F= 0
  DO i=101,200
    IF (i/=5 .AND. i/=6 .AND. i/=9) THEN
      !units 5-6-9 are commonly reserved for console I/O
      INQUIRE(UNIT=i,OPENED=lopen,IOSTAT=ios)
      IF (ios==0) THEN
        IF (.NOT. lOpen) THEN
          F=i
		  !WRITE(*,*) "Get Unit F=", F
          RETURN
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  WRITE(*,*) "Get Unit F=", F
  RETURN

ENDSUBROUTINE GetUnit

SUBROUTINE Str_Append(Str,Length,C)
!.when length of str <Length, append character C (default for C = "_")
  IMPLICIT NONE
  CHARACTER*(*),INTENT(INOUT)      :: Str !,S
  INTEGER,      INTENT(IN)         :: Length
  CHARACTER(1), INTENT(IN),OPTIONAL:: C
  DO
    IF(LEN_TRIM(Str)>=Length) EXIT
    IF(PRESENT(C)) THEN; Str=TRIM(Str)//C
    ELSE; Str=TRIM(Str)//"_"
    ENDIF
  ENDDO
  !IF(LEN_TRIM(Str)>Length)  ...!!!
  RETURN
ENDSUBROUTINE Str_Append

SUBROUTINE Str_Resize(Str,Length)
  IMPLICIT NONE
  CHARACTER*(*),INTENT(INOUT):: Str !,S
  INTEGER,      INTENT(IN)   :: Length
  INTEGER:: L
  !
  L=LEN_TRIM(Str)
  IF(LEN_TRIM(Str)>Length) Str(Length+1:L)=" "
  RETURN
ENDSUBROUTINE Str_Resize

INTEGER FUNCTION CarToInt(C) !conversion of char C to integer, IF(C NOT Valid) CarToInt=999
  IMPLICIT NONE
  CHARACTER,INTENT(IN):: C
  INTEGER  :: I
  CarToInt=999
  I=ICHAR(C)
  IF (I>=ICHAR('0').AND.I<=ICHAR('9')) CarToInt= I - ICHAR('0')
ENDFUNCTION CarToInt

!!! ! Internal read gives a variable string-represented numbers
!!! CHARACTER*12 str
!!! str = '123456'
!!! READ (str,'(i6)') i

SUBROUTINE WrdToInt(StrIn,IOut) !conversion of string StrIn to integer IOut
  IMPLICIT NONE
  CHARACTER*(*),INTENT(IN) :: StrIn
  INTEGER,      INTENT(OUT):: iOut
  !
  INTEGER:: Pos,LenStr, I, iSign
  !
  Pos=1; IOut=0; iSign=1; LenStr=LEN(StrIn)
  DO WHILE (POS<=LenStr.AND.StrIn(POS:POS)==' '); POS=POS+1; ENDDO !find first non blank
  DO WHILE (POS<=LenStr.AND.StrIn(Pos:Pos)/=' ') !do until find EoL or first blank
    IF (StrIn(Pos:Pos)=='-') THEN; ISIGN=-1
    ELSE
      I=ICHAR(StrIn(POS:POS))
      IF (I>=ICHAR('0').AND.I<=ICHAR('9')) THEN; IOut=  IOut*10 + I - ICHAR('0')
      ELSE;                                      IOut=0; RETURN
      ENDIF
    ENDIF
    Pos=Pos+1
  ENDDO
  IOut=iSIGN*IOut
  RETURN
ENDSUBROUTINE WrdToInt

SUBROUTINE WrdToInt2(StrIn,IOut,IOk)
!.conversion of string StrIn to integer IOut, with IO check
  IMPLICIT NONE
  CHARACTER*(*),INTENT(IN) :: StrIn
  INTEGER,      INTENT(OUT):: iOut
  !
  LOGICAL:: IOk
  INTEGER:: Pos,LenStr, I, iSign
  !
  Pos=1; IOut=0; iSign=1; LenStr=LEN(StrIn); IOk=.TRUE.
  DO WHILE (POS<=LenStr.AND.StrIn(POS:POS)==' '); POS=POS+1; ENDDO !find first non blank
  DO WHILE (POS<=LenStr.AND.StrIn(Pos:Pos)/=' ') !do until find EoL or first blank
    IF (StrIn(Pos:Pos)=='-') THEN; ISIGN=-1
    ELSE
      I=ICHAR(StrIn(POS:POS))
      IF (I>=ICHAR('0').AND.I<=ICHAR('9')) THEN; IOut=  IOut*10 + I - ICHAR('0')
      ELSE;                                      IOk=.FALSE.; RETURN
      ENDIF
    ENDIF
    Pos=Pos+1
  ENDDO
  IOut=iSIGN*IOut
  RETURN
ENDSUBROUTINE WrdToInt2

SUBROUTINE WrdToReal(StrIn,ROut) !conversion of string StrIn to real ROut
  IMPLICIT NONE
  CHARACTER*(*),INTENT(IN) :: StrIn
  REAL(dp),     INTENT(OUT):: ROut
  !
  INTEGER :: Pos, PosSep !PosSep:decimal separator
  REAL(dp):: Expon
  INTEGER :: LenStr, I, iSIGN
  !
  LenStr=LEN(StrIn)
  ROut=Zero; iSIGN=1; Expon=0; PosSep=0; Pos=1
  DO WHILE (Pos<=LenStr.AND.StrIn(Pos:Pos)==' '); Pos=Pos+1; ENDDO !find 1st non-blank
  DO WHILE ((StrIn(Pos:Pos)/=' ') & !___!First Read Mantissa Part
       .AND.(StrIn(Pos:Pos)/='D') &
       .AND.(StrIn(Pos:Pos)/='E') &
       .AND.(Pos<=LenStr))
    IF (StrIn(Pos:Pos)=='-')     THEN; iSIGN=-1
    ELSEIF (StrIn(Pos:Pos)=='.') THEN; PosSep=Pos
    ELSE
      I=ICHAR(StrIn(Pos:Pos))
      IF((I>=ICHAR('0')).AND.(I<=ICHAR('9'))) THEN
        ROut=  ROut*10 + I - ICHAR('0')
      ELSE
        ROut=Zero; RETURN
      ENDIF
    ENDIF
    Pos=Pos+1
  ENDDO
  ROut=iSIGN*ROut
  IF ((PosSep>0).AND.(PosSep<Pos-1)) THEN
    DO I=PosSep+1,Pos-1; ROut=ROut/10; ENDDO
  ENDIF
  iSIGN=1
  EXPON=0
  IF((Pos<=LenStr) &
  .AND.((StrIn(Pos:Pos)=='D').OR.(StrIn(Pos:Pos)=='d') &
    .OR.(StrIn(Pos:Pos)=='E').OR.(StrIn(Pos:Pos)=='e'))) THEN
    DO !WHILE ((StrIn(Pos:Pos)/=' ').AND.(Pos<=LenStr))
      Pos=Pos+1
      IF (Pos>LenStr)          EXIT
      IF (StrIn(Pos:Pos)==' ') EXIT
      IF (StrIn(Pos:Pos)=='-') THEN; iSIGN=-1; CYCLE; ENDIF
      IF (StrIn(Pos:Pos)=='+') THEN; iSIGN=+1; CYCLE; ENDIF
      I=ICHAR(StrIn(Pos:Pos))
      IF ((I>=ICHAR('0')).AND.(I<=ICHAR('9'))) THEN
        EXPON= EXPON*10 + I-ICHAR('0')
      ELSE
        ROut=Zero; RETURN
      ENDIF
    ENDDO
  ENDIF
  IF(EXPON>0) ROut=ROut*EXP(iSIGN*EXPON*LOG(10D0))
  RETURN
ENDSUBROUTINE WrdToReal

SUBROUTINE ReadRValsV(Line,NRead,vX)
!reads array vX from string Line; index of last element read is returned in NRead
  IMPLICIT NONE
  CHARACTER(LEN=*),        INTENT(IN) ::Line
  INTEGER,                 INTENT(OUT):: NRead
  REAL(dp),DIMENSION(dimV),INTENT(OUT):: vX
  !
  CHARACTER(512):: L,Word
  LOGICAL       :: EoL
  INTEGER       :: N
  !
  L=TRIM(Line)
  vX(1:dimV)=Zero
  DO N=1,dimV
    CALL LinToWrd(L,Word,EoL)
    CALL WrdToReal(Word,vX(N))
    IF(EoL) EXIT
  ENDDO
  NRead=N
ENDSUBROUTINE ReadRValsV

!~ SUBROUTINE LinToWrd_(Line,Word,Eol) !Reads Words From A String, return UPPERCASE word
  !~ IMPLICIT NONE
  !~ LOGICAL,INTENT(OUT):: EoL !TRUE= EndOfLine reached, or "!" found, when EoL, Line="!"
  !~ INTEGER:: ILow,Lenc,INext,I,J
  !~ CHARACTER(LEN=*) Line !=Input/Output= String=  Words Separated By Spaces or Tabs
  !~ CHARACTER(LEN=*) Word !=OUTPUT
  !~ !
  !~ Lenc=LEN_TRIM(Line)
  !~ Line=TRIM(Line)//'!'
  !~ EoL=.FALSE.

  !~ ILow=0 !______________________search first non blank/tab/!
  !~ DO
    !~ ILow=ILow+1
    !~ IF(Line(ILow:ILow)=='!') THEN
      !~ Word=  '!'; EoL=  .TRUE.
      !~ RETURN !end of line found ->exit
    !~ END IF
    !~ IF(Line(ILow:ILow)/=' ' .AND. Line(ILow:ILow)/=T_) EXIT
  !~ ENDDO

  !~ INext=ILow !__________________iLow= 1st non(blank,tab,!) character
  !~ DO !__________________________search for the last contiguous character that's not (blank/tab/!)
    !~ INext=INext+1
    !~ IF (Line(INext:INext)=='!' &
    !~ .OR.Line(INext:INext)==' ' &
    !~ .OR.Line(INext:INext)==T_) EXIT
  !~ ENDDO
  !~ Word=Line(ILow:INext-1)
  !~ !
  !~ DO I=1,LEN_TRIM(Word) !_______________!Output= UpperCase !!!
    !~ J=ICHAR(Word(I:I))
    !~ !IF (65<=J.AND.J<=90) Word(I:I)=CHAR(J+32) !LowerCase
    !~ IF (97<=J.AND.J<=122) Word(I:I)=CHAR(J-32) !UpperCase !!!
  !~ ENDDO

  !~ ILow=INext
  !~ DO !then, go to 1st occurrence of non (blank,tab) char, or to EoL
    !~ IF (Line(ILow:ILow)=='!') THEN
      !~ Line= '!'; EoL= .TRUE.
      !~ RETURN !__remaining Line is empty
    !~ END IF
    !~ IF(Line(ILow:ILow)/=' ' .AND. Line(ILow:ILow)/=T_) EXIT
    !~ ILow=ILow+1
  !~ ENDDO !_______________________iLow= 1st non(blank,tab) character
  !~ !
  !~ Line=Line(ILow:Lenc)
  !~ !
  !~ RETURN
!~ ENDSUBROUTINE LinToWrd_

SUBROUTINE LinToWrd(Line,Word,Eol,Cod)
!==
!== Reads Words From A String, with Case control
!==
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(INOUT)::Line !=a line of Words Separated By Spaces or Tabs
  CHARACTER(LEN=*),INTENT(OUT)  ::Word !=First Word of Line
  LOGICAL,         INTENT(OUT)  :: EoL !=EndOfLine reached, or "!" found, if EoL then Line="!"
  CHARACTER(LEN=2),INTENT(IN),OPTIONAL:: Cod !="UP","LO","NO"
  !
  INTEGER:: ILow,Lenc,INext,I,J
  !
  Lenc= LEN_TRIM(Line)
  Line= TRIM(Line)//'!'
  EoL=.FALSE.
  !
  ILow=0
  !==< search first non blank/tab/!
  DO
    ILow=ILow+1
    IF(Line(ILow:ILow)=='!') THEN
      Word=  '!'; EoL=  .TRUE.; RETURN !end of line found ->exit
    END IF
    IF(Line(ILow:ILow)/=' ' .AND. Line(ILow:ILow)/=T_) EXIT
  ENDDO
  !
  INext=ILow ! iLow= 1st non(blank,tab,!) character
  !==< search for the last contiguous character that's not (blank/tab/!)
  DO
    INext=INext+1
    IF (Line(INext:INext)=='!' &
    .OR.Line(INext:INext)==' ' &
    .OR.Line(INext:INext)==T_) EXIT
  ENDDO
  Word=Line(ILow:INext-1)
  !
  IF(PRESENT(Cod)) THEN
    !"LO"-> lower case
    !"UP"-> upper case
    !other value -> do not modify character case
    IF(Cod=="LO".OR.Cod=="UP") THEN
      DO I=1,LEN_TRIM(Word)
        J=ICHAR(Word(I:I))
        IF(Cod=="LO".AND.65<=J.AND.J<=90) Word(I:I)=CHAR(J+32) !LowerCase
        IF(Cod=="UP".AND.97<=J.AND.J<=122) Word(I:I)=CHAR(J-32) !UpperCase !!!
      ENDDO
    ENDIF
    !
  ELSE !
    !
    !Cod not given -> upper case
    DO I=1,LEN_TRIM(Word)
      J=ICHAR(Word(I:I))
      IF(97<=J.AND.J<=122) Word(I:I)=CHAR(J-32) !defaut=UpperCase !!!
    ENDDO
    !
  ENDIF
  !
  ILow=INext
  !== then, go to 1st occurrence of non (blank,tab) char, or to EoL
  DO
    IF (Line(ILow:ILow)=='!') THEN
      Line= '!'
      EoL= .TRUE.
      RETURN !== remaining Line is empty
    END IF
    IF(Line(ILow:ILow)/=' ' .AND. Line(ILow:ILow)/=T_) EXIT
    ILow=ILow+1
  ENDDO
  !== iLow= 1st non(blank,tab) character
  Line=Line(ILow:Lenc)
  !
  RETURN
ENDSUBROUTINE LinToWrd

SUBROUTINE AppendToEnd(Line,Word,Eol)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(INOUT):: Line
  CHARACTER(LEN=*),INTENT(INOUT):: Word !=
  LOGICAL,         INTENT(INOUT):: EoL  !=
  CHARACTER(LEN=255):: WW
  IF(.NOT. Eol) THEN
    IF(Word=="END") THEN !when END followed by keyword, append
      CALL LinToWrd(Line,WW,EoL); Word=TRIM(Word)//TRIM(WW)
    ENDIF
  ENDIF
ENDSUBROUTINE AppendToEnd

SUBROUTINE In_Car(S,C) !Screen/Kbd dialog
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: S
  CHARACTER,       INTENT(OUT):: C
  !
  CHARACTER(LEN=1)::Word
  INTEGER  :: ist,nch,J
  !
  WRITE(*,'(A)',ADVANCE='NO') TRIM(S)
  READ (*,'(A)',ADVANCE='NO',SIZE=nch,IOSTAT=ist) Word
  J=ICHAR(Word(1:1))
  !IF (65<=J.AND.J<=90) Word(I:I)=CHAR(J+32) !-> LowerCase
  IF (97<=J.AND.J<=122) Word(1:1)=CHAR(J-32) !-> UpperCase
  C=Word(1:1)
ENDSUBROUTINE In_Car

FUNCTION In_Str(S) !Screen/Kbd dialog
  IMPLICIT NONE
  CHARACTER(LEN=3)           ::In_Str
  CHARACTER(LEN=*),INTENT(IN):: S
  !
  CHARACTER(LEN=1)::Word
  INTEGER  :: ist,nch,J
  !
  nch=3
  WRITE(*,'(A)',ADVANCE='NO') TRIM(S) !; READ(*,*) StrIn !'(A)',ADVANCE='NO') StrIn
  READ(*, '(A)',ADVANCE='NO',SIZE=nch,IOSTAT=ist) Word
  !!
  J=ICHAR(Word(1:1))
  !IF (65<=J.AND.J<=90) Word(I:I)=CHAR(J+32) !-> LowerCase
  IF (97<=J.AND.J<=122) Word(1:1)=CHAR(J-32) !-> UpperCase
  !!
  in_Str=TRIM(Word)
ENDFUNCTION In_Str

FUNCTION In_Int(S,iMin,iMax) !Screen/Kbd dialog
  IMPLICIT NONE
  CHARACTER(LEN=*):: S
  INTEGER         :: In_Int,iMin,iMax
  INTEGER         :: ist,nch
  !
  CHARACTER(32):: StrIn
  INTEGER      :: iOut
  LOGICAL      :: IOk
  DO
    WRITE(*,'(A)',ADVANCE='NO') TRIM(S) !; READ(*,*) StrIn !'(A)',ADVANCE='NO') StrIn
    READ(*,'(A)',ADVANCE='NO',SIZE=nch,IOSTAT=ist) StrIn
    CALL WrdToInt2(StrIn,IOut,IOk) !iOk= StrIn is valid string
    IF(IOk .AND. IOut>=iMin .AND. IOut<=iMax) EXIT
  ENDDO
  in_Int=iOut
ENDFUNCTION In_Int

FUNCTION In_Real(S,rMin,rMax) !Screen/Kbd dialog
  !USE M_Glob;
  IMPLICIT NONE
  CHARACTER(LEN=32):: S
  REAL(dp)        :: In_Real, r,rMin,rMax
  DO
    WRITE(*,'(A)',ADVANCE='NO') TRIM(S); READ *,r
    IF(r>=rMin .AND. r<=rMax) EXIT
  ENDDO
  In_Real=r
ENDFUNCTION In_Real

SUBROUTINE IntToStr3(I,Str) !conversion I<10000 to string(4), for writing index
  IMPLICIT NONE
  INTEGER         ::I
  CHARACTER(LEN=3),INTENT(OUT)::Str
  INTEGER::J
  IF(I>=1000.OR.I<0) THEN; Str="999"; RETURN; ENDIF
  !-> INTERNAL WRITE STATEMENT
  WRITE(Str,'(I3)') I
  DO J=1,LEN(Str); IF(Str(J:J)==' ') Str(J:J)='0'; ENDDO
ENDSUBROUTINE IntToStr3

CHARACTER(LEN=3) FUNCTION FIntToStr3(I) !conversion I<10000 to string(4), for writing index
  IMPLICIT NONE
  INTEGER,INTENT(IN)::I
  !
  CHARACTER(LEN=3)::Str
  INTEGER::J
  IF(I>=1000.OR.I<0) THEN; Str="999"; RETURN; ENDIF
  !-> INTERNAL WRITE STATEMENT
  WRITE(Str,'(I3)') I
  DO J=1,LEN(Str); IF(Str(J:J)==' ') Str(J:J)='0'; ENDDO
  FIntToStr3=Str
ENDFUNCTION FIntToStr3

SUBROUTINE IntToStr4(I,Str) !conversion I<10000 to string(4), for writing index
  IMPLICIT NONE
  INTEGER,         INTENT(IN) ::I
  CHARACTER(LEN=4),INTENT(OUT)::Str
  INTEGER::J
  IF(I>=10000.OR.I<0) THEN; Str="9999"; RETURN; ENDIF
  WRITE(Str,'(I4)') I
  DO J=1,LEN(Str); IF(Str(J:J)==' ') Str(J:J)='0'; ENDDO
ENDSUBROUTINE IntToStr4

CHARACTER(LEN=4) FUNCTION FIntToStr4(I) !conversion I<10000 to string(4), for writing index
  IMPLICIT NONE
  INTEGER,INTENT(IN)::I
  !
  CHARACTER(LEN=4)::Str
  INTEGER::J
  IF(I>=10000.OR.I<0) THEN; Str="9999"; RETURN; ENDIF
  !-> INTERNAL WRITE STATEMENT
  WRITE(Str,'(I4)') I
  DO J=1,LEN(Str); IF(Str(J:J)==' ') Str(J:J)='0'; ENDDO
  FIntToStr4=Str
ENDFUNCTION FIntToStr4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ModIOWrite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ModIOWrite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ModIOWrite
SUBROUTINE OutStrVec( &
& fOut,Vec, & !in
& Opt_I,Opt_J,Opt_S,Opt_R,Opt_C,vL,CR) !in, optional
!.write Counter Opt_I, tab, string Str, tab, array Vec,
!.according to filter vL and format Opt_C="E"/"F"/"G"
  IMPLICIT NONE
  INTEGER,              INTENT(IN):: fOut
  REAL(dp),DIMENSION(:),INTENT(IN):: Vec
  INTEGER,              INTENT(IN),OPTIONAL:: Opt_I,Opt_J !counter(s)
  CHARACTER(LEN=*),     INTENT(IN),OPTIONAL:: Opt_S !a string
  REAL(dp),             INTENT(IN),OPTIONAL:: Opt_R !values
  CHARACTER,            INTENT(IN),OPTIONAL:: Opt_C !format
  LOGICAL, DIMENSION(:),INTENT(IN),OPTIONAL:: vL !write or not, should be same size as Vec !!!
  LOGICAL,              INTENT(IN),OPTIONAL:: CR
  !
  LOGICAL,DIMENSION(SIZE(Vec))::vLocal !write or not
  INTEGER:: J
  LOGICAL:: CR_
  !
  IF(PRESENT(CR)) THEN; CR_= CR
  ELSE                ; CR_= .TRUE.
  ENDIF
  IF(PRESENT(vL)) THEN; vLocal= vL
  ELSE                ; vLocal= .TRUE.
  ENDIF
  IF(PRESENT(Opt_I)) WRITE(fOut,'(I7,A1)',ADVANCE='NO') Opt_I,T_
  IF(PRESENT(Opt_J)) WRITE(fOut,'(I7,A1)',ADVANCE='NO') Opt_J,T_
  IF(PRESENT(Opt_S)) WRITE(fOut,'(A,A1)',ADVANCE='NO') TRIM(Opt_S),T_
  IF(PRESENT(Opt_R)) WRITE(fOut,'(G15.6,A1)',ADVANCE='NO') Opt_R,T_
  DO J=1,SIZE(Vec)
    IF(vLocal(J)) THEN
      IF(Vec(J)==Zero) THEN; WRITE(fOut,'(14X,A1,A1)',ADVANCE='NO') "0", T_
      ELSE
        IF(PRESENT(Opt_C)) THEN
        SELECT CASE(Opt_C)
          CASE("E"); WRITE(fOut,'( E15.6,A1)',ADVANCE='NO')  Vec(J), T_
          CASE("F"); WRITE(fOut,'( F15.6,A1)',ADVANCE='NO')  Vec(J), T_
          CASE("G"); WRITE(fOut,'( G15.6,A1)',ADVANCE='NO')  Vec(J), T_
          CASE("7"); WRITE(fOut,'( F7.1, A1)',ADVANCE='NO')  Vec(J), T_
        ENDSELECT
        ELSE
          WRITE(fOut,'( G15.6,A1)',ADVANCE='NO')  Vec(J), T_
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  IF(CR_) WRITE(fOut,*)
ENDSUBROUTINE OutStrVec

SUBROUTINE Str_Upper(Str)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(INOUT):: Str
  !
  INTEGER:: I,J,D
  !
  D= ICHAR('A') - ICHAR('a')
  DO I=1,LEN(Str)
    J= ICHAR(Str(I:I))
    IF(J>=ICHAR('a') .AND. J<=ICHAR('z')) Str(I:I)= CHAR(J+D)
  ENDDO
  RETURN
ENDSUBROUTINE Str_Upper

!~ CHARACTER(LEN=*) FUNCTION FStr_Upper(Str)
  !~ IMPLICIT NONE
  !~ CHARACTER(LEN=*), INTENT(IN):: Str
  !~ !
  !~ INTEGER:: I,J,D
  !~ CHARACTER(LEN=LEN(Str)):: UpStr
  !~ !
  !~ D= ICHAR('A') - ICHAR('a')
  !~ UpStr= TRIM(Str)
  !~ DO I=1,LEN(Str)
    !~ J= ICHAR(Str(I:I))
    !~ IF(J>=ICHAR('a') .AND. J<=ICHAR('z')) UpStr(I:I)= CHAR(J+D)
  !~ ENDDO
  !~ FStr_Upper= TRIM(UpStr)
  !~ !
  !~ RETURN
!~ END FUNCTION FStr_Upper

!---

  SUBROUTINE OpenFile(F, FILE)
   IMPLICIT NONE
   INTEGER :: F
   CHARACTER(LEN=*) :: FILE
   !write(*,*) "Opening File Unit F=", F, "=", trim(FILE)
   OPEN(F,FILE=FILE)
  END SUBROUTINE

   !---

  SUBROUTINE CloseFile(F)
    IMPLICIT NONE
    INTEGER:: F
    !write(*,*) "Closing File Unit F=", F
    CLOSE(F)
  END SUBROUTINE

ENDMODULE M_IoTools

