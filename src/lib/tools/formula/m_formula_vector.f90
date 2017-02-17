module M_Formula_Vector

  !------------------------------------------------------
  ! Fortran Translation of a Matlab Recursive Code 
  ! Formula Parser
  !------------------------------------------------------
  use m_trace,   only : iDebug
  use m_iotools, only : WrdtoReal
  use m_kinds,   only : dp
  implicit none
  private
  real(dp) :: ppcm_div
  public :: Formula_Vector_Read
  
contains 

  !---

  subroutine Formula_Vector_Read ( formula, name, s, b, isok, div )
    implicit none
    character(len=*), intent(in) :: formula
    character(len=*), intent(in) :: name(:) 
    logical, intent(out) :: isok   
    !---
    real(dp), intent(out) :: b(:)
    real(dp), intent(out), optional :: div 
    character(len=*), intent(out) :: s
    !---
    real(dp) :: factor = 1.D0
    b = 0.D0
    !--== Call Recursive Subroutine
    ppcm_div = 1.D0
    call  get_formula_vector ( formula, name, factor, s, b, isok )
    if (present(div)) div = ppcm_div
  
  end subroutine Formula_Vector_Read

  !---------------------------------------------------------------------

  recursive subroutine get_formula_vector ( formula, name, factor, s, b, isok )
    implicit none
    character(len=*), intent(in) :: formula
    character(len=*), intent(in) :: name(:) 
    real(dp) :: factor
    
    !---
    real(dp), intent(inout) :: b(:)
    character(len=*), intent(out) :: s
    logical, intent(out) :: isok   
    integer :: n
    integer :: End_s, i
    character(len=20)  :: coef, sym
    character(len=1)   :: next
    character(len=len(formula)) :: scopy 
    !---
    n = size(name)
    isok = .true.
    s = formula;
    coef  = ''                  !% initialize stoichiometric coefficient
    sym   = ''                  !% initialize chemical symbol
    !
    !// store the ppm of the denominator
    !
    do while(.not.(isempty(s)))                         !% continue while string is non-empty
      !
      End_s    = len(trim(s));
      next     = s(End_s:End_s);                       !% pick next token from end of string
      !
      s        = s(1:End_s-1);                         !% decrement the formula string
      if(iDebug==5) write(*,*) "Next String, ", Next, ' ', trim(s) !
      if (isletter(next)) then                         !% token is a letter a-z or A-Z
        sym= trim(next)//trim(sym)                     !% build a chemical symbol (He, Mg, etc.)
        call strmatch(sym,name, i)                     !% test against known symbols
        if (anyy(i)) then                              !% SYM is a recognized symbol
          if(iDebug==5) write(*,*) "Match Name :", name(i) !%
          b(i) = b(i) + factor*stoichiometry(coef);    !% increment stoich. coeff
          coef = '';                                   !% reset the stoichiometric coefficient
          sym  = '';                                   !% reset the chemical symbol
        end if                                         !%
      elseif (isnumber(next)) then                     !% token is a period or a number 0-9
        coef= trim(next)//trim(coef)                   !% build a stoichiometric coefficient
        if(iDebug==5) write(*,*) "Is number ", coef
        call emptystring(sym, isok);                   !% test that SYM is empty
      elseif (isright(next)) then                      !% token is a right paranthesis )
        if(iDebug==5) write(*,*) "Is Right ", next
        !%-----------------------------------------------
        !% call myself recursive + increase premultiplic. factor
        scopy = trim(s)
        call get_formula_vector(scopy,name, factor*stoichiometry(coef), s, b, isok); 
        if(iDebug==5) write(*,*) 
        if (.not.(isok)) exit
        !%-----------------------------------------------
        coef   = '';                                   !% reset the stoichiometric coefficient
        call emptystring(sym, isok);                   !% test that SYM is empty
      elseif (isleft(next))   then                     !% token a left paranthesis (
        if(iDebug==5) write(*,*) "Is Left ", next
        call emptystring(sym, isok);                   !% test that SYM is empty
        return;                                        !% return (from recursive call only)
      end if
    end do
    !
    call emptystring(sym, isok);                       !% test that SYM is empty
    ! 
  end subroutine get_formula_vector

  !---

  !% test that sym is empty
  subroutine emptystring(sym, isok)                      
    implicit none
    character(len=*) :: sym
    logical :: isok
    !write(*,*) 'Testing string: ',sym
    if (.not.(trim(sym).eq.'')) then
       !write(*,*) 'Unrecognized string: ',sym
       isok=.false. 
    else
       isok =.true.
    end if
  end subroutine emptystring

  !---

  ! %Verify that input string is a period or a number in the range 0-9 
  function isnumber(s) result(Ok)
    implicit none
    character(len=1) :: s
    logical :: Ok
    Ok = ( index( '0123456789./', s) > 0 )
  end function isnumber

  !%Verify that input string is a letter 
  function isletter(s) result(Ok)
    implicit none
    character(len=1) :: s
    logical :: Ok
    integer :: idx
    idx = index( 'abcdefghijklmnopqrstuvwxyz-+ABCDEFGHIJKLMNOPQRSTUVWXYZ', s) 
    Ok = (idx > 0)
  end function isletter

  !%Verify that input string is a left paranthesis (
  function isleft(s) result(Ok)
    implicit none
    character(len=1) :: s
    logical :: Ok
    Ok = s .eq. '('
  end function isleft

  !%Verify that input string is a left paranthesis (
  function isempty(s) result(Ok)
    implicit none
    character(len=1) :: s
    logical :: Ok
    Ok =  (trim(s).eq.'')
  end function isempty

  !%Verify that input string is a left paranthesis (
  function isright(s) result(Ok)
    implicit none
    character(len=1) :: s
    logical :: Ok
    Ok =  s.eq. ')'
  end function isright

  !%Convert input string to number (empty input string is treated as 1)
  function anyy(i) result(Ok)
    implicit none
    integer :: i
    logical :: Ok
    Ok = (.not.(i.eq.0))
  end function anyy

  !%Convert input string to number (empty input string is treated as 1)
  function stoichiometry(s) result(n)
    implicit none
    character(len=*) :: s
    real(dp) :: n 
    if (isempty(s)) then
       n = 1.d0;
    else
       call str2num(s, n);
    end if
  end function stoichiometry

  !%Convert input string to number (empty input string is treated as 1)
  subroutine str2num(s, x)
    ! convert fraction a/b to a real value
    implicit none
    character(len=*) :: s
    real(dp) :: x
    real(dp) :: y 
    call WrdToRealFraction(s,x, y)
    x = x / y;
    !write(*,*) "x,y,ppcm", x,y, ppcm_div
    ppcm_div = ppcm_div*y
  end subroutine str2num

  subroutine WrdToRealFraction(sinput, a, b)
    ! convert fraction a/b to a real value
    implicit none
    character(len=*) :: sinput
    character(len=60) :: s, sa, sb
    real(dp), intent(out) :: a,b 
    !--
    integer :: idx
    integer :: slen
    
    !--
    s=adjustl(trim(sinput))
    slen = len(s)
    !--
    idx = index(s, '/')
    
    if (idx.eq.0) then
       call WrdToReal(s,a)
       b = 1.D0
    else 
       if (idx.eq.1) then
          a = 1.D0
       else
          sa = s(1:idx-1)
          call WrdToReal(sa,a)
       end if
       
       if (idx.eq.slen) then
          b = 1.D0
       else
          sb = s(idx+1:slen)
          call WrdToReal(sb,b)
       end if
       
    end if
    !write(*,*) "a,b=", a, b
  
  end subroutine WrdToRealFraction

  
  !%Convert input string to number (empty input string is treated as 1)
  subroutine strmatch(s, name, idx) 
    implicit none
    integer :: idx
    character(len=*) :: s
    character(len=*) :: name(:)
    integer :: i, n
    !write(*,*) "Search Match For ", trim(s)
    n = size(name)
    idx = 0
    do i=1, n
       if (trim(name(i)).eq.trim(s)) then
          idx = i
          exit
       end if
    end do
    if (idx>0 .and. iDebug==5) write(*,*) "Matching ", trim(name(idx))

  end subroutine strmatch

end module M_Formula_Vector


!//================= MATLAB ORIGINAL FILE ====================================

!!$%Adds stoichiometric coefficients for a specified  chemical  species using 
!!$%recursive string parsing.Eg [s,b]=get_formula_vector('(CH3(CH2)2COOH)1.5'
!!$%,[0 0 0],{'C','H','O'},1)  returns s='' and b=[6 12 3]. The chemical ele-
!!$%ments need not be specified in any particular order, but their order must 
!!$%coincide with the stoichiometric coefficients also specified in the input
!!$%
!!$%Author         : Tore Haug-Warberg
!!$%Version        : March 10 2000 (THW)
!!$%
!!$%function [s,b] = get_formula_vector(s,b,name,factor)
!!$%
!!$%        s      = Chemical formula [string].
!!$%        b      = Stoichiometric coefficients [vector].
!!$%        name   = Recognized formula symbols [cell array].
!!$%        factor = External multiplier.
!!$%
!!$
!!$ function [s,b,isok] = get_formula_vector(sx,bx,name,factor)
!!$ 
!!$ isok = true;
!!$ s = sx;
!!$ b = bx;
!!$ coef       = '';                  % initialize stoichiometric coefficient
!!$ sym        = '';                             % initialize chemical symbol
!!$
!!$ while not(isempty(s))                % continue while string is non-empty
!!$   next     = s(end);                 % pick next token from end of string
!!$   s        = s(1:end-1);                   % decrement the formula string
!!$   if isletter(next)                        % token is a letter a-z or A-Z
!!$     sym    = strcat(next,sym);   % build a chemical symbol (He, Mg, etc.)
!!$     i      = strmatch(sym,name,'exact');     % test against known symbols
!!$     if anyy(i)                                % SYM is a recognized symbol
!!$       disp(char(name(i)))
!!$       b(i) = b(i) + factor*stoichiometry(coef); % increment stoich. coeff
!!$       coef = '';                   % reset the stoichiometric coefficient
!!$       sym  = '';                              % reset the chemical symbol
!!$     end                                                                 %
!!$   elseif isnumber(next)               % token is a period or a number 0-9
!!$     coef   = strcat(next,coef);      % build a stoichiometric coefficient
!!$     [sym,isok]    = emptystring(sym);                   % test that SYM is empty
!!$   elseif isright(next)                   % token is a right paranthesis )
!!$     [s,b]  = get_formula_vector(s,b,name, ...   % call myself recursively
!!$              factor*stoichiometry(coef)); % increase premultiplic. factor
!!$     coef   = '';                   % reset the stoichiometric coefficient
!!$     [sym,isok]    = emptystring(sym);                   % test that SYM is empty
!!$   elseif isleft(next)                        % token a left paranthesis (
!!$     [sym,isok]    = emptystring(sym);                   % test that SYM is empty
!!$     return;                           % return (from recursive call only)
!!$   end                                                                   %
!!$ end                                                                     %
!!$
!!$ [sym,isok] = emptystring(sym);                          % test that SYM is empty
!!$     
!!$%Verify that input string is empty
!!$ function [s,isok] = emptystring(s)
!!$ if (not(isempty(s)))
!!$   disp(['Unrecognized string: ',s]); 
!!$   isok = false
!!$ else
!!$   isok = true
!!$ end
!!$ sym = '';
!!$ 
!!$%Verify that input string is a period or a number in the range 0-9
!!$ function [i] = isnumber(s)
!!$ i = or(and(ge(abs(s),48),le(abs(s),57)),eq(abs(s),46));
!!$
!!$%Verify that input string is a left paranthesis (
!!$ function [i] = isleft(s)
!!$ i = eq(abs(s),40);
!!$
!!$%Verify that input string is a right paranthesis )
!!$ function [i] = isright(s)
!!$ i = eq(abs(s),41);
!!$ 
!!$%Convert input string to number (empty input string is treated as 1)
!!$ function [n] = stoichiometry(s)
!!$ if isempty(s)
!!$   n = 1;
!!$ else
!!$   n = str2num(s);
!!$ end

!//================= end MATLAB ORIGINAL FILE ====================================
