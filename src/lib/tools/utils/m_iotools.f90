!.module for basic output procedures
module M_IOTools
  use M_Kinds
  implicit none
  !
  private
  !
  public:: cUpper, cLower
  public:: Str_Append
  public:: Str_Upper
  public:: CarToInt,WrdToInt
  public:: WrdToReal,ReadRValsV
  public:: IntToStr3,IntToStr4,FIntToStr3,FIntToStr4
  public:: In_Int,In_Car,In_Str
  public:: LinToWrd
  public:: AppendToEnd
  public:: OutStrVec
  public:: GetUnit
  public:: OpenFile
  public:: CloseFile
  !
  integer,parameter,public:: dimV=32 !max.dimension for routine ReadRValsV
  !
  character,parameter:: T_= char(9)

contains

character function cUpper(c)
  character,intent(in):: c
  integer:: j,d
  d= ichar('A') - ichar('a')
  cUpper= c
  j= ichar(c)
  if(j>=ichar('a') .and. j<=ichar('z')) cUpper= char(j+d)
end function cUpper

character function cLower(c)
  character,intent(in):: c
  integer:: j,d
  d= ichar('A') - ichar('a')
  cLower= c
  j= ichar(c)
  if(j>=ichar('A') .and. j<=ichar('Z')) cLower= char(j-d)
end function cLower

subroutine GetUnit(F) !returns a free unit number.
!-- Author: John Burkardt
!-- A "free" FORTRAN unit number is an integer between 1 and 99
!-- which is not currently associated with an I/O device.
!-- A free FORTRAN unit number is needed in order to open a file 
!-- with the open command.
!-- IUNIT=0
!-- no free FORTRAN unit could be found, although all 99 units were checked
  implicit none
  integer,intent(out)::F
  integer:: i, ios
  logical:: lOpen
  F= 0
  do i=10,200
    !if (i/=5 .and. i/=6 .and. i/=9) then
    !units 5-6-9 are commonly reserved for console I/O
    inquire(unit=i,opened=lopen,iostat=ios)
    if (ios==0) then
      if (.not. lOpen) then
        F=i
        return
      end if
    end if
    !end if
  end do
  return
end subroutine GetUnit

subroutine Str_Append(Str,Length,C)
!.when length of str <Length, append character C (default for C = "_")
  implicit none
  character*(*),intent(inout)      :: Str !,S
  integer,      intent(in)         :: Length
  character(1), intent(in),optional:: C
  do
    if(len_trim(Str)>=Length) exit
    if(present(C)) then; Str=trim(Str)//C
    else; Str=trim(Str)//"_"
    end if
  end do
  !if(len_trim(Str)>Length)  ...!!!
  return
end subroutine Str_Append

subroutine Str_Resize(Str,Length)
  implicit none
  character*(*),intent(inout):: Str !,S
  integer,      intent(in)   :: Length
  integer:: L
  !
  L=len_trim(Str)
  if(len_trim(Str)>Length) Str(Length+1:L)=" "
  return
end subroutine Str_Resize

integer function CarToInt(C)
!conversion of char C to integer, if(C NOT Valid) CarToInt=999
  implicit none
  character,intent(in):: C
  integer  :: I
  CarToInt=999
  I=ichar(C)
  if (I>=ichar('0').and.I<=ichar('9')) CarToInt= I - ichar('0')
end function CarToInt

!!! ! Internal read gives a variable string-represented numbers
!!! character*12 str
!!! str = '123456'
!!! read (str,'(i6)') i

subroutine WrdToInt(StrIn,IOut) !conversion of string StrIn to integer IOut
  implicit none
  character*(*),intent(in) :: StrIn
  integer,      intent(out):: iOut
  !
  integer:: Pos,LenStr, I, iSign
  !
  Pos=1; IOut=0; iSign=1; LenStr=len(StrIn)
  do while (POS<=LenStr.and.StrIn(POS:POS)==' ')
    POS=POS+1
  end do !find first non blank
  do while (POS<=LenStr.and.StrIn(Pos:Pos)/=' ') !do until find EoL or first blank
    if (StrIn(Pos:Pos)=='-') then
      ISIGN=-1
    else
      I=ichar(StrIn(POS:POS))
      if (I>=ichar('0').and.I<=ichar('9')) then; IOut=  IOut*10 + I - ichar('0')
      else;                                      IOut=0; return
      end if
    end if
    Pos=Pos+1
  end do
  IOut=iSIGN*IOut
  return
end subroutine WrdToInt

subroutine WrdToInt2(StrIn,    IOut,IOk)
!.conversion of string StrIn to integer IOut, with IO check
  implicit none
  character*(*),intent(in) :: StrIn
  integer,      intent(out):: iOut
  logical,      intent(out):: IOk
  !
  integer:: Pos,LenStr, I, iSign
  !
  IOk=.true.
  Pos=1; IOut=0; iSign=1; LenStr=len(StrIn)
  do while (POS<=LenStr.and.StrIn(POS:POS)==' ')
  !-------------------------------------------------find first non blank
    POS=POS+1
  end do 
  do while (POS<=LenStr.and.StrIn(Pos:Pos)/=' ')
  !-------------------------------------do until find EoL or first blank
    if (StrIn(Pos:Pos)=='-') then; ISIGN=-1
    else
      I=ichar(StrIn(POS:POS))
      if (I>=ichar('0').and.I<=ichar('9')) then
        IOut=  IOut*10 + I - ichar('0')
      else
        IOk=.false.
        return
      end if
    end if
    Pos=Pos+1
  end do
  IOut=iSIGN*IOut
  return
end subroutine WrdToInt2

subroutine WrdToReal(StrIn,    ROut)
!-- conversion of string StrIn to real ROut
  implicit none
  character*(*),intent(in) :: StrIn
  real(dp),     intent(out):: ROut
  !
  integer :: Pos, PosSep !PosSep:decimal separator
  real(dp):: Expon
  integer :: LenStr, I, iSIGN
  !
  LenStr=len(StrIn)
  ROut=Zero; iSIGN=1; Expon=0; PosSep=0; Pos=1
  do while (Pos<=LenStr.and.StrIn(Pos:Pos)==' '); Pos=Pos+1; end do !find 1st non-blank
  do while ((StrIn(Pos:Pos)/=' ') & !___!First Read Mantissa Part
       .and.(StrIn(Pos:Pos)/='D') &
       .and.(StrIn(Pos:Pos)/='E') &
       .and.(Pos<=LenStr))
    if (StrIn(Pos:Pos)=='-')     then; iSIGN=-1
    elseif (StrIn(Pos:Pos)=='.') then; PosSep=Pos
    else
      I=ichar(StrIn(Pos:Pos))
      if((I>=ichar('0')).and.(I<=ichar('9'))) then
        ROut=  ROut*10 + I - ichar('0')
      else
        ROut=Zero; return
      end if
    end if
    Pos=Pos+1
  end do
  ROut=iSIGN*ROut
  if ((PosSep>0).and.(PosSep<Pos-1)) then
    do I=PosSep+1,Pos-1; ROut=ROut/10; end do
  end if
  iSIGN=1
  EXPON=0
  if((Pos<=LenStr) &
  .and.((StrIn(Pos:Pos)=='D').or.(StrIn(Pos:Pos)=='d') &
    .or.(StrIn(Pos:Pos)=='E').or.(StrIn(Pos:Pos)=='e'))) then
    do !while ((StrIn(Pos:Pos)/=' ').and.(Pos<=LenStr))
      Pos=Pos+1
      if (Pos>LenStr)          exit
      if (StrIn(Pos:Pos)==' ') exit
      if (StrIn(Pos:Pos)=='-') then; iSIGN=-1; cycle; end if
      if (StrIn(Pos:Pos)=='+') then; iSIGN=+1; cycle; end if
      I=ichar(StrIn(Pos:Pos))
      if ((I>=ichar('0')).and.(I<=ichar('9'))) then
        EXPON= EXPON*10 + I-ichar('0')
      else
        ROut=Zero; return
      end if
    end do
  end if
  if(EXPON>0) ROut=ROut*exp(iSIGN*EXPON*log(10D0))
  return
end subroutine WrdToReal

subroutine ReadRValsV(Line,NRead,vX)
!-- reads array vX from string Line; 
!-- index of last element read is returned in NRead
  implicit none
  character(len=*),        intent(in) ::Line
  integer,                 intent(out):: NRead
  real(dp),dimension(dimV),intent(out):: vX
  !
  character(512):: L,Word
  logical       :: EoL
  integer       :: N
  !
  L=trim(Line)
  vX(1:dimV)=Zero
  do N=1,dimV
    call LinToWrd(L,Word,EoL)
    call WrdToReal(Word,vX(N))
    if(EoL) exit
  end do
  NRead=N
end subroutine ReadRValsV

!! subroutine LinToWrd_(Line,Word,Eol) !Reads Words From A String, return UPPERcase word
!!   implicit none
!!   logical,intent(out):: EoL !TRUE= EndOfLine reached, or "!" found, when EoL, Line="!"
!!   integer:: ILow,Lenc,INext,I,J
!!   character(len=*) Line !=Input/Output= String=  Words Separated By Spaces or Tabs
!!   character(len=*) Word !=OUTPUT
!!   !
!!   Lenc=len_trim(Line)
!!   Line=trim(Line)//'!'
!!   EoL=.false.
!! 
!!   ILow=0 !______________________search first non blank/tab/!
!!   do
!!     ILow=ILow+1
!!     if(Line(ILow:ILow)=='!') then
!!       Word=  '!'; EoL=  .true.
!!       return !end of line found ->exit
!!     end if
!!     if(Line(ILow:ILow)/=' ' .and. Line(ILow:ILow)/=T_) exit
!!   end do
!! 
!!   INext=ILow !__________________iLow= 1st non(blank,tab,!) character
!!   do !__________________________search for the last contiguous character that's not (blank/tab/!)
!!     INext=INext+1
!!     if (Line(INext:INext)=='!' &
!!     .or.Line(INext:INext)==' ' &
!!     .or.Line(INext:INext)==T_) exit
!!   end do
!!   Word=Line(ILow:INext-1)
!!   !
!!   do I=1,len_trim(Word) !_______________!Output= UpperCase !!!
!!     J=ichar(Word(I:I))
!!     !if (65<=J.and.J<=90) Word(I:I)=char(J+32) !LowerCase
!!     if (97<=J.and.J<=122) Word(I:I)=char(J-32) !UpperCase !!!
!!   end do
!! 
!!   ILow=INext
!!   do !then, go to 1st occurrence of non (blank,tab) char, or to EoL
!!     if (Line(ILow:ILow)=='!') then
!!       Line= '!'; EoL= .true.
!!       return !__remaining Line is empty
!!     end if
!!     if(Line(ILow:ILow)/=' ' .and. Line(ILow:ILow)/=T_) exit
!!     ILow=ILow+1
!!   end do !_______________________iLow= 1st non(blank,tab) character
!!   !
!!   Line=Line(ILow:Lenc)
!!   !
!!   return
!! end subroutine LinToWrd_

subroutine LinToWrd(Line,Word,Eol,Cod)
!--
!-- Reads Words From A String, with Case control
!--
  implicit none
  character(len=*),intent(inout)::Line !=a line of Words Separated By Spaces or Tabs
  character(len=*),intent(out)  ::Word !=First Word of Line
  logical,         intent(out)  :: EoL !=EndOfLine reached, or "!" found, if EoL then Line="!"
  character(len=2),intent(in),optional:: Cod !="UP","LO","NO"
  !
  integer:: ILow,Lenc,INext,I,J
  !
  Lenc= len_trim(Line)
  Line= trim(Line)//'!'
  EoL=.false.
  !
  ILow=0
  !----------------------------------------search first non blank or tab
  do
    ILow=ILow+1
    if(Line(ILow:ILow)=='!') then
      Word=  '!'; EoL=  .true.; return !end of line found ->exit
    end if
    if(Line(ILow:ILow)/=' ' .and. Line(ILow:ILow)/=T_) exit
  end do
  !--------------------------------------------------------------------/
  !
  INext=ILow ! iLow= 1st non(blank,tab,!) character
  !
  !-----search for the last contiguous character that's not blank or tab
  do
    INext=INext+1
    if (Line(INext:INext)=='!' &
    .or.Line(INext:INext)==' ' &
    .or.Line(INext:INext)==T_) exit
  end do
  Word=Line(ILow:INext-1)
  !--------------------------------------------------------------------/
  !
  if(present(Cod)) then
    !"LO"-> lower case
    !"UP"-> upper case
    !other value -> do not modify character case
    if(Cod=="LO".or.Cod=="UP") then
      do I=1,len_trim(Word)
        J=ichar(Word(I:I))
        if(Cod=="LO".and.65<=J.and.J<=90) Word(I:I)=char(J+32) !LowerCase
        if(Cod=="UP".and.97<=J.and.J<=122) Word(I:I)=char(J-32) !UpperCase !!!
      end do
    end if
    !
  else !
    !
    !Cod not given -> upper case
    do I=1,len_trim(Word)
      J=ichar(Word(I:I))
      if(97<=J.and.J<=122) Word(I:I)=char(J-32) !defaut=UpperCase !!!
    end do
    !
  end if
  !
  ILow=INext
  !--------then, go to 1st occurrence of non (blank,tab) char, or to EoL
  do
    if (Line(ILow:ILow)=='!') then
      Line= '!'
      EoL= .true.
      return !-----------------------------------remaining Line is empty
    end if
    if(Line(ILow:ILow)/=' ' .and. Line(ILow:ILow)/=T_) exit
    ILow=ILow+1
  end do
  !-----------------------------------iLow= 1st non(blank,tab) character
  Line=Line(ILow:Lenc)
  !
  return
end subroutine LinToWrd

subroutine AppendToEnd(Line,Word,Eol)
  implicit none
  character(len=*),intent(inout):: Line
  character(len=*),intent(inout):: Word !=
  logical,         intent(inout):: EoL  !=
  character(len=255):: WW
  if(.not. Eol) then
    if(Word=="END") then !when end followed by keyword, append
      call LinToWrd(Line,WW,EoL)
      Word=trim(Word)//trim(WW)
    end if
  end if
end subroutine AppendToEnd

subroutine In_Car(S,C) !Screen/Kbd dialog
  implicit none
  character(len=*),intent(in) :: S
  character,       intent(out):: C
  !
  character(len=1)::Word
  integer  :: ist,nch,J
  !
  write(*,'(A)',advance="no") trim(S)
  read (*,'(A)',advance="no",size=nch,iostat=ist) Word
  J=ichar(Word(1:1))
  !if (65<=J.and.J<=90) Word(I:I)=char(J+32) !-> LowerCase
  if (97<=J.and.J<=122) Word(1:1)=char(J-32) !-> UpperCase
  C=Word(1:1)
end subroutine In_Car

function In_Str(S) !Screen/Kbd dialog
  implicit none
  character(len=3)           ::In_Str
  character(len=*),intent(in):: S
  !
  character(len=1)::Word
  integer  :: ist,nch,J
  !
  nch=3
  write(*,'(A)',advance="no") trim(S) !; read(*,*) StrIn !'(A)',advance="no") StrIn
  read(*, '(A)',advance="no",size=nch,iostat=ist) Word
  !!
  J=ichar(Word(1:1))
  !if (65<=J.and.J<=90) Word(I:I)=char(J+32) !-> LowerCase
  if (97<=J.and.J<=122) Word(1:1)=char(J-32) !-> UpperCase
  !!
  in_Str=trim(Word)
end function In_Str

function In_Int(S,iMin,iMax) !Screen/Kbd dialog
  implicit none
  character(len=*):: S
  integer         :: In_Int,iMin,iMax
  integer         :: ist,nch
  !
  character(32):: StrIn
  integer      :: iOut
  logical      :: IOk
  do
    write(*,'(A)',advance="no") trim(S) !; read(*,*) StrIn !'(A)',advance="no") StrIn
    read(*,'(A)',advance="no",size=nch,iostat=ist) StrIn
    call WrdToInt2(StrIn,IOut,IOk) !iOk= StrIn is valid string
    if(IOk .and. IOut>=iMin .and. IOut<=iMax) exit
  end do
  in_Int=iOut
end function In_Int

function In_Real(S,rMin,rMax) !Screen/Kbd dialog
  !use M_Glob;
  implicit none
  character(len=32):: S
  real(dp)        :: In_Real, r,rMin,rMax
  do
    write(*,'(A)',advance="no") trim(S); read *,r
    if(r>=rMin .and. r<=rMax) exit
  end do
  In_Real=r
end function In_Real

subroutine IntToStr3(I,Str) !conversion I<10000 to string(4), for writing index
  implicit none
  integer         ::I
  character(len=3),intent(out)::Str
  integer::J
  if(I>=1000.or.I<0) then; Str="999"; return; end if
  !-> INTERNAL write STATEMENT
  write(Str,'(I3)') I
  do J=1,len(Str); if(Str(J:J)==' ') Str(J:J)='0'; end do
end subroutine IntToStr3

character(len=3) function FIntToStr3(I) !conversion I<10000 to string(4), for writing index
  implicit none
  integer,intent(in)::I
  !
  character(len=3)::Str
  integer::J
  if(I>=1000.or.I<0) then; Str="999"; return; end if
  !-> INTERNAL write STATEMENT
  write(Str,'(I3)') I
  do J=1,len(Str); if(Str(J:J)==' ') Str(J:J)='0'; end do
  FIntToStr3=Str
end function FIntToStr3

subroutine IntToStr4(I,Str) !conversion I<10000 to string(4), for writing index
  implicit none
  integer,         intent(in) ::I
  character(len=4),intent(out)::Str
  integer::J
  if(I>=10000.or.I<0) then; Str="9999"; return; end if
  write(Str,'(I4)') I
  do J=1,len(Str); if(Str(J:J)==' ') Str(J:J)='0'; end do
end subroutine IntToStr4

character(len=4) function FIntToStr4(I) !conversion I<10000 to string(4), for writing index
  implicit none
  integer,intent(in)::I
  !
  character(len=4)::Str
  integer::J
  if(I>=10000.or.I<0) then; Str="9999"; return; end if
  !-> INTERNAL write STATEMENT
  write(Str,'(I4)') I
  do J=1,len(Str); if(Str(J:J)==' ') Str(J:J)='0'; end do
  FIntToStr4=Str
end function FIntToStr4

subroutine OutStrVec( &
& fOut,Vec, & !in
& Opt_I,Opt_J,Opt_S,Opt_R,Opt_C,vL,CR) !in, optional
!-- write Counter Opt_I, tab, string Str, tab, array Vec,
!-- according to filter vL and format Opt_C="E"/"F"/"G"
  !
  implicit none
  !
  integer,              intent(in):: fOut
  real(dp),dimension(:),intent(in):: Vec
  integer,              intent(in),optional:: Opt_I,Opt_J !counter(s)
  character(len=*),     intent(in),optional:: Opt_S !a string
  real(dp),             intent(in),optional:: Opt_R !values
  character,            intent(in),optional:: Opt_C !format
  logical, dimension(:),intent(in),optional:: vL !write or not, should be same size as Vec !!!
  logical,              intent(in),optional:: CR
  !
  logical,dimension(size(Vec))::vLocal !write or not
  integer:: J
  logical:: CR_
  !
  if(present(CR)) then  ; CR_= CR
  else                  ; CR_= .true.
  end if
  !
  if(present(vL)) then  ; vLocal= vL
  else                  ; vLocal= .true.
  end if
  !
  if(present(Opt_I)) write(fOut,'(I7,A1)',advance="no") Opt_I,T_
  if(present(Opt_J)) write(fOut,'(I7,A1)',advance="no") Opt_J,T_
  if(present(Opt_S)) write(fOut,'(A,A1)',advance="no") trim(Opt_S),T_
  if(present(Opt_R)) write(fOut,'(G15.6,A1)',advance="no") Opt_R,T_
  do J=1,size(Vec)
    if(vLocal(J)) then
      if(Vec(J)==Zero) then; write(fOut,'(14X,A1,A1)',advance="no") "0", T_
      else
        if(present(Opt_C)) then
        select case(Opt_C)
          case("E"); write(fOut,'( E15.6,A1)',advance="no")  Vec(J), T_
          case("F"); write(fOut,'( F15.6,A1)',advance="no")  Vec(J), T_
          case("G"); write(fOut,'( G15.6,A1)',advance="no")  Vec(J), T_
          case("7"); write(fOut,'( F7.1, A1)',advance="no")  Vec(J), T_
        end select
        else
          write(fOut,'( G15.6,A1)',advance="no")  Vec(J), T_
        end if
      end if
    end if
  end do
  if(CR_) write(fOut,*)
end subroutine OutStrVec

subroutine Str_Upper(Str)
  implicit none
  character(len=*), intent(inout):: Str
  !
  integer:: I,J,D
  !
  D= ichar('A') - ichar('a')
  do I=1,len(Str)
    J= ichar(Str(I:I))
    if(J>=ichar('a') .and. J<=ichar('z')) Str(I:I)= char(J+D)
  end do
  return
end subroutine Str_Upper

!! character(len=*) function FStr_Upper(Str)
  !! implicit none
  !! character(len=*), intent(in):: Str
  !! !
  !! integer:: I,J,D
  !! character(len=len(Str)):: UpStr
  !! !
  !! D= ichar('A') - ichar('a')
  !! UpStr= trim(Str)
  !! do I=1,len(Str)
    !! J= ichar(Str(I:I))
    !! if(J>=ichar('a') .and. J<=ichar('z')) UpStr(I:I)= char(J+D)
  !! end do
  !! FStr_Upper= trim(UpStr)
  !! !
  !! return
!! end function FStr_Upper

!---

  subroutine OpenFile(F, FILE)
   implicit none
   integer :: F
   character(len=*) :: FILE
   !write(*,*) "Opening File Unit F=", F, "=", trim(FILE)
   open(F,file=FILE)
  end subroutine

   !---

  subroutine CloseFile(F)
    implicit none
    integer:: F
    !write(*,*) "Closing File Unit F=", F
    close(F)
  end subroutine

end module M_IoTools

