MODULE M_VarStr
!.module for string input, extracted from ISO_VARYING_STRING
!
!ISO_VARYING_STRING=
!Written by J.L.Schonfelder 
!Incorporating suggestions by C.Tanasescu, C.Weber, J.Wagener and W.Walter,
!and corrections due to L.Moss, M.Cohen, P.GrIFiths, B.T.Smith
!and many other members of the committee ISO/IEC JTC1/SC22/WG5

PRIVATE 

PUBLIC:: T_VarStr, VarStr_Get, VarStr_Char, S_Ass_C !, LEN
  
TYPE T_VarStr
 PRIVATE 
 CHARACTER,DIMENSION(:),POINTER :: chars 
ENDTYPE T_VarStr 
  
CHARACTER,PARAMETER :: blank= " " 
  
!----- GENERIC PROCEDURE INTERFACE DEFINITIONS -------------------------------! 
  
!----- LEN interface ---------------------------------------------------------! 
INTERFACE LEN 
  MODULE PROCEDURE len_s   ! length of string
ENDINTERFACE 
  
!----- Conversion procedure interfaces ---------------------------------------!
INTERFACE VAR_STR
  MODULE PROCEDURE c_to_s   ! character to string
ENDINTERFACE 
  
INTERFACE VarStr_Char
  MODULE PROCEDURE &
  & s_to_c, &   ! string to character
  & s_to_fix_c  ! string to specified length character
ENDINTERFACE 
  
!----- ASSIGNMENT interfaces -------------------------------------------------! 
INTERFACE ASSIGNMENT(=) 
  MODULE PROCEDURE &
  & s_ass_s, &   ! string = string
  & c_ass_s, &   ! character= string
  & s_ass_c      ! string = character
ENDINTERFACE 
  
!----- Concatenation operator interfaces -------------------------------------!
INTERFACE OPERATOR(//) 
  MODULE PROCEDURE &
  & s_concat_s, &  ! string//string
  & s_concat_c, &  ! string//character
  & c_concat_s     ! character//string
ENDINTERFACE 
  
!----- Repeated Concatenation interface --------------------------------------! 
INTERFACE REPEAT 
  MODULE PROCEDURE repeat_s
ENDINTERFACE 

!----- Input procedure interfaces --------------------------------------------!
INTERFACE VarStr_GET
  MODULE PROCEDURE &
  & get_d_eor !, &    ! default unit, EoR termination
  !& get_u_eor,    & ! specified unit, EoR termination
  !& get_d_tset_s, & ! default unit, string set termination
  !& get_u_tset_s, & ! specified unit, string set termination
  !& get_d_tset_c, & ! default unit, char set termination
  !& get_u_tset_c    ! specified unit, char set termination
ENDINTERFACE 
  
CONTAINS

!----- LEN Procedure ---------------------------------------------------------! 
FUNCTION len_s(string) 
  TYPE(T_VarStr),INTENT(IN) :: string 
  INTEGER                         :: len_s 
  ! returns the length of the string argument or zero iFin there iStr no current 
  ! string value 
  IF(.NOT.ASSOCIATED(string%chars))THEN; len_s= 0 
  ELSE;                                  len_s= SIZE(string%chars) 
  ENDIF 
ENDFUNCTION len_s 
  
!----- Conversion Procedures ------------------------------------------------! 
FUNCTION c_to_s(chr) 
  TYPE(T_VarStr)              :: c_to_s 
  CHARACTER(LEN=*),INTENT(IN) :: chr 
  ! returns the string consisting of the characters char 
  INTEGER                     :: lc,i 
  lc=LEN(chr) 
  ALLOCATE(c_to_s%chars(1:lc)) 
  DO i=1,lc; c_to_s%chars(i)= chr(i:i); ENDDO 
ENDFUNCTION c_to_s 
  
FUNCTION s_to_c(string) 
! returns the characters of string as an automatically sized character 
  TYPE(T_VarStr),INTENT(IN)   :: string 
  CHARACTER(LEN=SIZE(string%chars)) :: s_to_c 
  INTEGER:: lc, i 
  lc=SIZE(string%chars) 
  DO i=1,lc ; s_to_c(i:i)= string%chars(i); ENDDO 
ENDFUNCTION s_to_c 
  
FUNCTION s_to_fix_c(string,length)
! returns the character of fixed length, length,
! containing the characters of string
! either padded with blanks or truncated on the right to fit
  TYPE(T_VarStr),INTENT(IN) :: string
  INTEGER,INTENT(IN)              :: length
  CHARACTER(LEN=length)           :: s_to_fix_c
  INTEGER                         :: lc, i
  lc=MIN(SIZE(string%chars),length)
  DO i=1,lc ; s_to_fix_c(i:i)= string%chars(i); ENDDO 
  IF(lc < length) & ! result longer than string padding needed
  & s_to_fix_c(lc+1:length)= blank
ENDFUNCTION s_to_fix_c

!----- ASSIGNMENT Procedures -------------------------------------------------!
SUBROUTINE s_ass_s(var,expr) 
  TYPE(T_VarStr),INTENT(OUT) :: var 
  TYPE(T_VarStr),INTENT(IN)  :: expr 
  !  assign a string value to a string variable overriding default assignement 
  !  reallocates string variable to size of string value and copies characters 
  ALLOCATE(var%chars(1:LEN(expr)))
  var%chars= expr%chars 
ENDSUBROUTINE s_ass_s 
 
SUBROUTINE c_ass_s(var,expr) 
  CHARACTER(LEN=*),INTENT(OUT)    :: var 
  TYPE(T_VarStr),INTENT(IN) :: expr 
  ! assign a string value to a character variable 
  ! iFin the string iStr longer than the character truncate the string on the right 
  ! iFin the string iStr shorter the character iStr blank padded on the right 
  INTEGER                         :: lc,ls,i 
  lc= LEN(var); ls= MIN(LEN(expr),lc) 
  DO i= 1,ls;    var(i:i)= expr%chars(i); ENDDO 
  DO i= ls+1,lc; var(i:i)= blank        ; ENDDO 
ENDSUBROUTINE c_ass_s 
  
SUBROUTINE s_ass_c(var,expr)
  TYPE(T_VarStr),INTENT(OUT) :: var 
  CHARACTER(LEN=*),INTENT(IN):: expr 
  !  assign a character value to a string variable 
  !  disassociates the string variable from its current value, allocates new 
  !  space to hold the characters and copies them from the character value 
  !  into this space.
  INTEGER                          :: lc, i 
  lc= LEN(expr) 
  ALLOCATE(var%chars(1:lc)) 
  DO i= 1,lc; var%chars(i)= expr(i:i); ENDDO 
ENDSUBROUTINE s_ass_c 
  
!----- Concatenation operator procedures ------------------------------------! 
FUNCTION s_concat_s(string_a,string_b)  ! string//string 
  TYPE(T_VarStr),INTENT(IN) :: string_a,string_b 
  TYPE(T_VarStr)            :: s_concat_s 
  INTEGER                         :: la,lb
  la= LEN(string_a); lb= LEN(string_b)
  ALLOCATE(s_concat_s%chars(1:la+lb)) 
  s_concat_s%chars(1:la)= string_a%chars 
  s_concat_s%chars(1+la:la+lb)= string_b%chars 
ENDFUNCTION s_concat_s 
  
FUNCTION s_concat_c(string_a,string_b)  ! string//character
  TYPE(T_VarStr),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  TYPE(T_VarStr)            :: s_concat_c 
  INTEGER                         :: la,lb, i 
  la= LEN(string_a); lb= LEN(string_b)
  ALLOCATE(s_concat_c%chars(1:la+lb)) 
  s_concat_c%chars(1:la)= string_a%chars 
  DO i= 1,lb; s_concat_c%chars(la+i)= string_b(i:i); ENDDO 
ENDFUNCTION s_concat_c 
  
FUNCTION c_concat_s(string_a,string_b)  ! character//string 
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  TYPE(T_VarStr),INTENT(IN) :: string_b 
  TYPE(T_VarStr)            :: c_concat_s 
  INTEGER                         :: la,lb, i 
  la= LEN(string_a); lb= LEN(string_b)
  ALLOCATE(c_concat_s%chars(1:la+lb)) 
  DO i= 1,la; c_concat_s%chars(i)= string_a(i:i); ENDDO 
  c_concat_s%chars(1+la:la+lb)= string_b%chars 
ENDFUNCTION c_concat_s 
  
!----- Reapeated concatenation procedures -----------------------------------! 
FUNCTION repeat_s(string,ncopies)                                     
 TYPE(T_VarStr),INTENT(IN) :: string 
 INTEGER,INTENT(IN)              :: ncopies 
 TYPE(T_VarStr)            :: repeat_s
 ! Returns a string produced by the concatenation of ncopies of the 
 ! argument string 
 INTEGER                         :: lr,ls, i 
 IF(ncopies < 0) THEN 
     WRITE(*,*) " Negative ncopies requested in REPEAT" 
     STOP 
 ENDIF 
 ls= LEN(string); lr= ls*ncopies 
 ALLOCATE(repeat_s%chars(1:lr))
 DO i= 1,ncopies 
   repeat_s%chars(1+(i-1)*ls:i*ls)= string%chars 
 ENDDO 
ENDFUNCTION repeat_s 
  
!----- Input string procedure -----------------------------------------------!
SUBROUTINE get_d_eor(string,maxlen,iostat) !,message)
! reads string from the default unit starting at next character in the file
! and terminating at the end of record or after maxlen characters.
  TYPE(T_VarStr),INTENT(OUT) :: string
  ! the string variable 
  ! to be filled with characters read from the file connected to the default unit
  INTEGER,INTENT(IN),OPTIONAL      :: maxlen
  ! if present, indicates the maximum number of characters that will be read from the file
  INTEGER,INTENT(OUT),OPTIONAL     :: iostat
  ! if present, used to return the status of the data transfer
  ! if absent errors cause termination
  !TYPE(T_VarStr),INTENT(IN),OPTIONAL:: message
  !
  CHARACTER(LEN=80) :: buffer
  INTEGER           :: ist,nch,toread,nb
  !! IF(PRESENT(message)) 
  IF(PRESENT(maxlen))THEN; toread=maxlen
  ELSE;                    toread=HUGE(1)
  ENDIF
  string= "" ! clears return string
  DO
    !repeatedly read buffer and add to string until EoR or maxlen reached
    IF(toread <= 0)EXIT
    nb=MIN(80,toread)
    READ(*,FMT='(A)',ADVANCE='NO',EOR=9999,SIZE=nch,IOSTAT=ist) buffer(1:nb)
    IF( ist /= 0 )THEN 
      IF(PRESENT(iostat)) THEN
        iostat=ist 
        RETURN 
      ELSE 
        WRITE(*,*) " Error No.",ist, &
        &          " during READ_STRING of varying string on default unit"
        STOP 
      ENDIF 
    ENDIF 
    string= string //buffer(1:nb)
    toread= toread - nb
  ENDDO
  IF(PRESENT(iostat)) iostat= 0
  RETURN
  9999 string= string //buffer(1:nch) 
  IF(PRESENT(iostat)) iostat= ist
ENDSUBROUTINE get_d_eor
  
ENDMODULE M_VarStr

