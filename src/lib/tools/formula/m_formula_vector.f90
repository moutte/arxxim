MODULE M_Formula_Vector

  !------------------------------------------------------
  ! Fortran Translation of a Matlab Recursive Code 
  ! Formula Parser
  !------------------------------------------------------
  USE m_trace,   ONLY : iDebug
  USE m_iotools, ONLY : WrdtoReal
  IMPLICIT NONE
  PRIVATE
  REAL(kind=8) :: ppcm_div
  PUBLIC :: Formula_Vector_Read
  
CONTAINS 

  !---

  SUBROUTINE Formula_Vector_Read ( formula, name, s, b, isok, div )
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: formula
    CHARACTER(LEN=*), INTENT(IN) :: name(:) 
    LOGICAL, INTENT(OUT) :: isok   
    !---
    REAL(kind=8), INTENT(OUT) :: b(:)
    REAL(kind=8), INTENT(OUT), optional :: div 
    CHARACTER(LEN=*), INTENT(OUT) :: s
    !---
    REAL(kind=8) :: factor = 1
    b = 0.D0
    !==== Call Recursive Subroutine
    ppcm_div = 1.D0
    CALL  get_formula_vector ( formula, name, factor, s, b, isok )
    IF (present(div)) div = ppcm_div
  
  END SUBROUTINE Formula_Vector_Read

  !===============================================================================================

  recursive SUBROUTINE get_formula_vector ( formula, name, factor, s, b, isok )
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: formula
    CHARACTER(LEN=*), INTENT(IN) :: name(:) 
    REAL(kind=8) :: factor
    
    !---
    REAL(kind=8), INTENT(INOUT) :: b(:)
    CHARACTER(LEN=*), INTENT(OUT) :: s
    LOGICAL, INTENT(OUT) :: isok   
    INTEGER :: n
    INTEGER :: End_s, i
    CHARACTER(LEN=20)  :: coef, sym
    CHARACTER(LEN=1)   :: next
    CHARACTER(LEN=LEN(formula)) :: scopy 
    !---
    n = SIZE(name)
    isok = .true.;
    s = formula;
    coef  = '';                                         !% initialize stoichiometric coefficient
    sym   = '';                                         !% initialize chemical symbol
    !
    !// store the ppm of the denominator
    !
    DO WHILE(.not.(isempty(s)))                         !% continue WHILE string is non-empty
      !
      End_s    = LEN(TRIM(s));
      next     = s(End_s:End_s);                       !% pick next token from END of string
      !
      s        = s(1:End_s-1);                         !% decrement the formula string
      IF(iDebug==5) WRITE(*,*) "Next String, ", Next, ' ', TRIM(s) !
      IF (isletter(next)) THEN                         !% token is a letter a-z or A-Z
        sym= TRIM(next)//TRIM(sym)                     !% build a chemical symbol (He, Mg, etc.)
        CALL strmatch(sym,name, i)                     !% test against known symbols
        IF (anyy(i)) THEN                              !% SYM is a recognized symbol
          IF(iDebug==5) WRITE(*,*) "Match Name :", name(i) !%
          b(i) = b(i) + factor*stoichiometry(coef);    !% increment stoich. coeff
          coef = '';                                   !% reset the stoichiometric coefficient
          sym  = '';                                   !% reset the chemical symbol
        END IF                                         !%
      ELSEIF (isnumber(next)) THEN                     !% token is a period or a number 0-9
        coef= TRIM(next)//TRIM(coef)                   !% build a stoichiometric coefficient
        IF(iDebug==5) WRITE(*,*) "Is number ", coef
        CALL emptystring(sym, isok);                   !% test that SYM is empty
      ELSEIF (isright(next)) THEN                      !% token is a right paranthesis )
        IF(iDebug==5) WRITE(*,*) "Is Right ", next
        !%-----------------------------------------------
        !% CALL myself recursive + increase premultiplic. factor
        scopy = TRIM(s)
        CALL get_formula_vector(scopy,name, factor*stoichiometry(coef), s, b, isok); 
        IF(iDebug==5) WRITE(*,*) 
        IF (.not.(isok)) EXIT
        !%-----------------------------------------------
        coef   = '';                                   !% reset the stoichiometric coefficient
        CALL emptystring(sym, isok);                   !% test that SYM is empty
      ELSEIF (isleft(next))   THEN                     !% token a left paranthesis (
        IF(iDebug==5) WRITE(*,*) "Is Left ", next
        CALL emptystring(sym, isok);                   !% test that SYM is empty
        RETURN;                                        !% RETURN (from recursive CALL ONLY)
      END IF
    END DO
    !
    CALL emptystring(sym, isok);                       !% test that SYM is empty
    ! 
  END SUBROUTINE get_formula_vector

  !---

  !% test that sym is empty
  SUBROUTINE emptystring(sym, isok)                      
    IMPLICIT NONE
    CHARACTER(LEN=*) :: sym
    LOGICAL :: isok
    !WRITE(*,*) 'Testing string: ',sym
    IF (.not.(TRIM(sym).eq.'')) THEN
       !WRITE(*,*) 'Unrecognized string: ',sym
       isok=.false. 
    ELSE
       isok =.true.
    END IF
  END SUBROUTINE emptystring

  !---

  ! %Verify that input string is a period or a number in the range 0-9 
  FUNCTION isnumber(s) result(Ok)
    IMPLICIT NONE
    CHARACTER(LEN=1) :: s
    LOGICAL :: Ok
    Ok = ( index( '0123456789./', s) > 0 )
  END FUNCTION isnumber

  !%Verify that input string is a letter 
  FUNCTION isletter(s) result(Ok)
    IMPLICIT NONE
    CHARACTER(LEN=1) :: s
    LOGICAL :: Ok
    INTEGER :: idx
    idx = index( 'abcdefghijklmnopqrstuvwxyz-+ABCDEFGHIJKLMNOPQRSTUVWXYZ', s) 
    Ok = (idx > 0)
  END FUNCTION isletter

  !%Verify that input string is a left paranthesis (
  FUNCTION isleft(s) result(Ok)
    IMPLICIT NONE
    CHARACTER(LEN=1) :: s
    LOGICAL :: Ok
    Ok = s .eq. '('
  END FUNCTION isleft

  !%Verify that input string is a left paranthesis (
  FUNCTION isempty(s) result(Ok)
    IMPLICIT NONE
    CHARACTER(LEN=1) :: s
    LOGICAL :: Ok
    Ok =  (TRIM(s).eq.'')
  END FUNCTION isempty

  !%Verify that input string is a left paranthesis (
  FUNCTION isright(s) result(Ok)
    IMPLICIT NONE
    CHARACTER(LEN=1) :: s
    LOGICAL :: Ok
    Ok =  s.eq. ')'
  END FUNCTION isright

  !%Convert input string to number (empty input string is treated as 1)
  FUNCTION anyy(i) result(Ok)
    IMPLICIT NONE
    INTEGER :: i
    LOGICAL :: Ok
    Ok = (.not.(i.eq.0))
  END FUNCTION anyy

  !%Convert input string to number (empty input string is treated as 1)
  FUNCTION stoichiometry(s) result(n)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: s
    REAL(kind=8) :: n 
    IF (isempty(s)) THEN
       n = 1.d0;
    ELSE
       CALL str2num(s, n);
    END IF
  END FUNCTION stoichiometry

  !%Convert input string to number (empty input string is treated as 1)
  SUBROUTINE str2num(s, x)
    ! convert fraction a/b to a REAL value
    IMPLICIT NONE
    CHARACTER(LEN=*) :: s
    REAL(kind=8) :: x
    REAL(kind=8) :: y 
    CALL WrdToRealFraction(s,x, y)
    x = x / y;
    !WRITE(*,*) "x,y,ppcm", x,y, ppcm_div
    ppcm_div = ppcm_div*y
  END SUBROUTINE str2num

  SUBROUTINE WrdToRealFraction(sinput, a, b)
    ! convert fraction a/b to a REAL value
    IMPLICIT NONE
    CHARACTER(LEN=*) :: sinput
    CHARACTER(LEN=60) :: s, sa, sb
    REAL(kind=8), INTENT(OUT) :: a,b 
    !--
    INTEGER :: idx
    INTEGER :: sLEN
    
    !--
    s=adjustl(TRIM(sinput))
    sLEN = LEN(s)
    !--
    idx = index(s, '/')
    
    IF (idx.eq.0) THEN
       CALL WrdToReal(s,a)
       b = 1.D0
    ELSE 
       IF (idx.eq.1) THEN
          a = 1.D0
       ELSE
          sa = s(1:idx-1)
          CALL WrdToReal(sa,a)
       END IF
       
       IF (idx.eq.sLEN) THEN
          b = 1.D0
       ELSE
          sb = s(idx+1:sLEN)
          CALL WrdToReal(sb,b)
       END IF
       
    END IF
    !WRITE(*,*) "a,b=", a, b
  
  END SUBROUTINE WrdToRealFraction

  
  !%Convert input string to number (empty input string is treated as 1)
  SUBROUTINE strmatch(s, name, idx) 
    IMPLICIT NONE
    INTEGER :: idx
    CHARACTER(LEN=*) :: s
    CHARACTER(LEN=*) :: name(:)
    INTEGER :: i, n
    !WRITE(*,*) "Search Match For ", TRIM(s)
    n = SIZE(name)
    idx = 0
    DO i=1, n
       IF (TRIM(name(i)).eq.TRIM(s)) THEN
          idx = i
          EXIT
       END IF
    END DO
    IF (idx>0 .and. iDebug==5) WRITE(*,*) "Matching ", TRIM(name(idx))

  END SUBROUTINE strmatch

ENDMODULE M_Formula_Vector


!//================= MATLAB ORIGINAL FILE ====================================

!!$%Adds stoichiometric coefficients for a specified  chemical  species using 
!!$%recursive string parsing.Eg [s,b]=get_formula_vector('(CH3(CH2)2COOH)1.5'
!!$%,[0 0 0],{'C','H','O'},1)  RETURNs s='' and b=[6 12 3]. The chemical ele-
!!$%ments need not be specified in any particular order, but their order must 
!!$%coincide with the stoichiometric coefficients also specified in the input
!!$%
!!$%Author         : Tore Haug-Warberg
!!$%Version        : March 10 2000 (THW)
!!$%
!!$%FUNCTION [s,b] = get_formula_vector(s,b,name,factor)
!!$%
!!$%        s      = Chemical formula [string].
!!$%        b      = Stoichiometric coefficients [vector].
!!$%        name   = Recognized formula symbols [cell array].
!!$%        factor = External multiplier.
!!$%
!!$
!!$ FUNCTION [s,b,isok] = get_formula_vector(sx,bx,name,factor)
!!$ 
!!$ isok = true;
!!$ s = sx;
!!$ b = bx;
!!$ coef       = '';                  % initialize stoichiometric coefficient
!!$ sym        = '';                             % initialize chemical symbol
!!$
!!$ WHILE not(isempty(s))                % continue while string is non-empty
!!$   next     = s(END);                 % pick next token from end of string
!!$   s        = s(1:END-1);                   % decrement the formula string
!!$   IF isletter(next)                        % token is a letter a-z or A-Z
!!$     sym    = strcat(next,sym);   % build a chemical symbol (He, Mg, etc.)
!!$     i      = strmatch(sym,name,'exact');     % test against known symbols
!!$     IF anyy(i)                                % SYM is a recognized symbol
!!$       disp(char(name(i)))
!!$       b(i) = b(i) + factor*stoichiometry(coef); % increment stoich. coeff
!!$       coef = '';                   % reset the stoichiometric coefficient
!!$       sym  = '';                              % reset the chemical symbol
!!$     END                                                                 %
!!$   ELSEIF isnumber(next)               % token is a period or a number 0-9
!!$     coef   = strcat(next,coef);      % build a stoichiometric coefficient
!!$     [sym,isok]    = emptystring(sym);                   % test that SYM is empty
!!$   ELSEIF isright(next)                   % token is a right paranthesis )
!!$     [s,b]  = get_formula_vector(s,b,name, ...   % CALL myself recursively
!!$              factor*stoichiometry(coef)); % increase premultiplic. factor
!!$     coef   = '';                   % reset the stoichiometric coefficient
!!$     [sym,isok]    = emptystring(sym);                   % test that SYM is empty
!!$   ELSEIF isleft(next)                        % token a left paranthesis (
!!$     [sym,isok]    = emptystring(sym);                   % test that SYM is empty
!!$     RETURN;                           % return (from recursive call only)
!!$   END                                                                   %
!!$ END                                                                     %
!!$
!!$ [sym,isok] = emptystring(sym);                          % test that SYM is empty
!!$     
!!$%Verify that input string is empty
!!$ FUNCTION [s,isok] = emptystring(s)
!!$ IF (not(isempty(s)))
!!$   disp(['Unrecognized string: ',s]); 
!!$   isok = false
!!$ ELSE
!!$   isok = true
!!$ END
!!$ sym = '';
!!$ 
!!$%Verify that input string is a period or a number in the range 0-9
!!$ FUNCTION [i] = isnumber(s)
!!$ i = or(and(ge(abs(s),48),le(abs(s),57)),eq(abs(s),46));
!!$
!!$%Verify that input string is a left paranthesis (
!!$ FUNCTION [i] = isleft(s)
!!$ i = eq(abs(s),40);
!!$
!!$%Verify that input string is a right paranthesis )
!!$ FUNCTION [i] = isright(s)
!!$ i = eq(abs(s),41);
!!$ 
!!$%Convert input string to number (empty input string is treated as 1)
!!$ FUNCTION [n] = stoichiometry(s)
!!$ IF isempty(s)
!!$   n = 1;
!!$ ELSE
!!$   n = str2num(s);
!!$ END

!//================= END MATLAB ORIGINAL FILE ====================================
