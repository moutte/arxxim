module M_VarStr
!.module for string input, extracted from ISO_VARYING_STRING
!
!ISO_VARYING_STRING=
!Written by J.L.Schonfelder 
!Incorporating suggestions by C.Tanasescu, C.Weber, J.Wagener and W.Walter,
!and corrections due to L.Moss, M.Cohen, P.Grifiths, B.T.Smith
!and many other members of the committee ISO/IEC JTC1/SC22/WG5

private 

public:: T_VarStr, VarStr_Get, VarStr_Char, S_Ass_C !, len
  
type T_VarStr
 private 
 character,dimension(:),pointer :: chars 
end type T_VarStr 
  
character,parameter :: blank= " " 
  
!----- generic procedure interface DEFINITIONS -------------------------------! 
  
!----- len interface ---------------------------------------------------------! 
interface len 
  module procedure len_s   ! length of string
endinterface 
  
!----- Conversion procedure interfaces ---------------------------------------!
interface VAR_STR
  module procedure c_to_s   ! character to string
endinterface 
  
interface VarStr_Char
  module procedure &
  & s_to_c, &   ! string to character
  & s_to_fix_c  ! string to specified length character
endinterface 
  
!----- assignMENT interfaces -------------------------------------------------! 
interface assignMENT(=) 
  module procedure &
  & s_ass_s, &   ! string = string
  & c_ass_s, &   ! character= string
  & s_ass_c      ! string = character
endinterface 
  
!----- Concatenation operator interfaces -------------------------------------!
interface operator(//) 
  module procedure &
  & s_concat_s, &  ! string//string
  & s_concat_c, &  ! string//character
  & c_concat_s     ! character//string
endinterface 
  
!----- Repeated Concatenation interface --------------------------------------! 
interface REPEAT 
  module procedure repeat_s
endinterface 

!----- Input procedure interfaces --------------------------------------------!
interface VarStr_GET
  module procedure &
  & get_d_eor !, &    ! default unit, EoR termination
  !& get_u_eor,    & ! specified unit, EoR termination
  !& get_d_tset_s, & ! default unit, string set termination
  !& get_u_tset_s, & ! specified unit, string set termination
  !& get_d_tset_c, & ! default unit, char set termination
  !& get_u_tset_c    ! specified unit, char set termination
endinterface 
  
contains

!----- len Procedure ---------------------------------------------------------! 
function len_s(string) 
  type(T_VarStr),intent(in) :: string 
  integer                         :: len_s 
  ! returns the length of the string argument or zero iFin there iStr no current 
  ! string value 
  if(.not.associateD(string%chars))then; len_s= 0 
  else;                                  len_s= size(string%chars) 
  end if 
end function len_s 
  
!----- Conversion Procedures ------------------------------------------------! 
function c_to_s(chr) 
  type(T_VarStr)              :: c_to_s 
  character(len=*),intent(in) :: chr 
  ! returns the string consisting of the characters char 
  integer                     :: lc,i 
  lc=len(chr) 
  allocate(c_to_s%chars(1:lc)) 
  do i=1,lc; c_to_s%chars(i)= chr(i:i); end do 
end function c_to_s 
  
function s_to_c(string) 
! returns the characters of string as an automatically sized character 
  type(T_VarStr),intent(in)   :: string 
  character(len=size(string%chars)) :: s_to_c 
  integer:: lc, i 
  lc=size(string%chars) 
  do i=1,lc ; s_to_c(i:i)= string%chars(i); end do 
end function s_to_c 
  
function s_to_fix_c(string,length)
! returns the character of fixed length, length,
! containing the characters of string
! either padded with blanks or truncated on the right to fit
  type(T_VarStr),intent(in) :: string
  integer,intent(in)              :: length
  character(len=length)           :: s_to_fix_c
  integer                         :: lc, i
  lc=MIN(size(string%chars),length)
  do i=1,lc ; s_to_fix_c(i:i)= string%chars(i); end do 
  if(lc < length) & ! result longer than string padding needed
  & s_to_fix_c(lc+1:length)= blank
end function s_to_fix_c

!----- assignMENT Procedures -------------------------------------------------!
subroutine s_ass_s(var,expr) 
  type(T_VarStr),intent(out) :: var 
  type(T_VarStr),intent(in)  :: expr 
  !  assign a string value to a string variable overriding default assignement 
  !  reallocates string variable to size of string value and copies characters 
  allocate(var%chars(1:len(expr)))
  var%chars= expr%chars 
end subroutine s_ass_s 
 
subroutine c_ass_s(var,expr) 
  character(len=*),intent(out)    :: var 
  type(T_VarStr),intent(in) :: expr 
  ! assign a string value to a character variable 
  ! iFin the string iStr longer than the character truncate the string on the right 
  ! iFin the string iStr shorter the character iStr blank padded on the right 
  integer                         :: lc,ls,i 
  lc= len(var); ls= min(len(expr),lc) 
  do i= 1,ls;    var(i:i)= expr%chars(i); end do 
  do i= ls+1,lc; var(i:i)= blank        ; end do 
end subroutine c_ass_s 
  
subroutine s_ass_c(var,expr)
  type(T_VarStr),intent(out) :: var 
  character(len=*),intent(in):: expr 
  !  assign a character value to a string variable 
  !  disassociates the string variable from its current value, allocates new 
  !  space to hold the characters and copies them from the character value 
  !  into this space.
  integer                          :: lc, i 
  lc= len(expr) 
  allocate(var%chars(1:lc)) 
  do i= 1,lc; var%chars(i)= expr(i:i); end do 
end subroutine s_ass_c 
  
!----- Concatenation operator procedures ------------------------------------! 
function s_concat_s(string_a,string_b)  ! string//string 
  type(T_VarStr),intent(in) :: string_a,string_b 
  type(T_VarStr)            :: s_concat_s 
  integer                         :: la,lb
  la= len(string_a); lb= len(string_b)
  allocate(s_concat_s%chars(1:la+lb)) 
  s_concat_s%chars(1:la)= string_a%chars 
  s_concat_s%chars(1+la:la+lb)= string_b%chars 
end function s_concat_s 
  
function s_concat_c(string_a,string_b)  ! string//character
  type(T_VarStr),intent(in) :: string_a 
  character(len=*),intent(in)     :: string_b 
  type(T_VarStr)            :: s_concat_c 
  integer                         :: la,lb, i 
  la= len(string_a); lb= len(string_b)
  allocate(s_concat_c%chars(1:la+lb)) 
  s_concat_c%chars(1:la)= string_a%chars 
  do i= 1,lb; s_concat_c%chars(la+i)= string_b(i:i); end do 
end function s_concat_c 
  
function c_concat_s(string_a,string_b)  ! character//string 
  character(len=*),intent(in)     :: string_a 
  type(T_VarStr),intent(in) :: string_b 
  type(T_VarStr)            :: c_concat_s 
  integer                         :: la,lb, i 
  la= len(string_a); lb= len(string_b)
  allocate(c_concat_s%chars(1:la+lb)) 
  do i= 1,la; c_concat_s%chars(i)= string_a(i:i); end do 
  c_concat_s%chars(1+la:la+lb)= string_b%chars 
end function c_concat_s 
  
!----- Reapeated concatenation procedures -----------------------------------! 
function repeat_s(string,ncopies)                                     
 type(T_VarStr),intent(in) :: string 
 integer,intent(in)              :: ncopies 
 type(T_VarStr)            :: repeat_s
 ! Returns a string produced by the concatenation of ncopies of the 
 ! argument string 
 integer                         :: lr,ls, i 
 if(ncopies < 0) then 
     write(*,*) " Negative ncopies requested in REPEAT" 
     stop 
 end if 
 ls= len(string); lr= ls*ncopies 
 allocate(repeat_s%chars(1:lr))
 do i= 1,ncopies 
   repeat_s%chars(1+(i-1)*ls:i*ls)= string%chars 
 end do 
end function repeat_s 
  
!----- Input string procedure -----------------------------------------------!
subroutine get_d_eor(string,maxlen,iostat) !,message)
! reads string from the default unit starting at next character in the file
! and terminating at the end of record or after maxlen characters.
  type(T_VarStr),intent(out) :: string
  ! the string variable 
  ! to be filled with characters read from the file connected to the default unit
  integer,intent(in),optional      :: maxlen
  ! if present, indicates the maximum number of characters that will be read from the file
  integer,intent(out),optional     :: iostat
  ! if present, used to return the status of the data transfer
  ! if absent errors cause termination
  !type(T_VarStr),intent(in),optional:: message
  !
  character(len=80) :: buffer
  integer           :: ist,nch,toread,nb
  !! if(present(message)) 
  if(present(maxlen))then; toread=maxlen
  else;                    toread=HUGE(1)
  end if
  string= "" ! clears return string
  do
    !repeatedly read buffer and add to string until EoR or maxlen reached
    if(toread <= 0)exit
    nb=MIN(80,toread)
    read(*,FMT='(A)',advance="no",EOR=9999,size=nch,iostat=ist) buffer(1:nb)
    if( ist /= 0 )then 
      if(present(iostat)) then
        iostat=ist 
        return 
      else 
        write(*,*) " Error No.",ist, &
        &          " during read_STRING of varying string on default unit"
        stop 
      end if 
    end if 
    string= string //buffer(1:nb)
    toread= toread - nb
  end do
  if(present(iostat)) iostat= 0
  return
  9999 string= string //buffer(1:nch) 
  if(present(iostat)) iostat= ist
end subroutine get_d_eor
  
end module M_VarStr

