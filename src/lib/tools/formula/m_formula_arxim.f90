module M_Formula_Arxim
!===============================================================
! Purpose :  Arxim Formulas Utils 
!---------------------------------------------------------------
! Arxim Formula Syntax 
!---------------------------------------------------------------
! CaCO3        = Ca(3)C(1)O(3) 
! HCO3-        = H(1)C(1)O(3)-(1)
! Ca+2         = Ca(1)+(2)
! CO3(MgCa)0.5 = C(2)O(6)MG(1)CA(1)/(2) 
!             or C(2)O(6)MG(1)CA(1)DIV(2)
!===============================================================
  use m_kinds
  implicit none

  private

  !// public functions
  public :: Formula_Arxim_Build
  public :: Formula_Arxim_Read
  
  interface Formula_Arxim_Build
    module procedure Formula_Arxim_Build_Real
    module procedure Formula_Arxim_Build_Int
  end interface

  logical, parameter :: DebFormula = .false.

contains

  subroutine Formula_Arxim_Build_Int ( vName, vCoef, OutFormula )
  !=====================================================================
  ! Purpose : Create an Arxim Formula from a table of Integer Coefficients
  ! Remark. Name can contain the pseudo-elements (+) and (-) 
  !=====================================================================
    implicit none
    character(len=*), intent(in) :: vName(:) 
    character(len=*), intent(out) :: OutFormula
    !---
    integer, intent(in) :: vCoef(:)
    integer :: n, i
    character(len=20) :: scoef
    !---
    OutFormula = ''
    n = size(vName)
    do i =1,n
      if (vCoef(i)==0) then 
      !// nothing to do
      else
        scoef=''
        write(scoef,'(I6)') vCoef(i)
        scoef = adjustl(scoef)
        OutFormula = trim(OutFormula)//trim(vName(i))//'('//trim(scoef)//')'
      end if
    end do
    !
  end subroutine Formula_Arxim_Build_Int

  !---
  
  subroutine Formula_Arxim_Build_Real ( vName, vCoef, OutFormula )
  !=====================================================================
  ! Purpose : Create an Arxim Formula from a table of Real Coefficients$
  ! Remark. Name can contain the pseudo-elements (+) and (-) 
  !=====================================================================
    character(len=*), intent(in) :: vName(:) 
    real(dp), intent(in) :: vCoef(:)
    character(len=*), intent(out) :: OutFormula
    !---
    integer :: n, i
    character(len=20) :: scoef
    real(dp) :: epsilon = 1.d-20
    !---
    OutFormula = ''
    n = size(vName)
    do i =1,n
      if (ABS(vCoef(i))<epsilon) then 
      !// nothing to do
      else
        !// write(*,*) trim(name(i)), ':', vCoef(i)
        scoef=''
        write(scoef,'(F10.6)') vCoef(i)
        scoef = adjustl(scoef)
        OutFormula = trim(OutFormula)//trim(vName(i))//'('//trim(scoef)//')'
      end if
    end do

  end subroutine Formula_Arxim_Build_Real

 
  !---

  subroutine Formula_Arxim_Read(sFormul,vNameEle,Z,nDiv,fOk,vStoik) 
  !=====================================================================
  !.read a chemical formula sFormul to a stoikiometric vector (vStoik)/Div 
  !.and the charge Z according to the element list vEle
  !.
  !.CAVEAT
  !.vStoik is an array of integers; nDiv the formula divisor
  !.thus, the "real" stoichio coeff for element vEle(i) is vstoik(i)/Div !!!
  !.
  !.-> used to read the elemental decompositon of a substance from its formula
  !=====================================================================
    use M_Formula_Utils
    use M_IOTools,only: WrdToInt
    !
    character(len=*),intent(in) :: sFormul       !the stoichio'formula of a substance
    character(len=*),intent(in) :: vNameEle(:)   !the current element name list
    integer,         intent(out):: Z,nDiv        !the charge and the "divisor"
    logical,         intent(out):: fOk           !the formula is consistent (with elemental basis)
    integer,         intent(out):: vStoik(:)     !the "formula vector"
    !
    character(len=80):: sRest
    character(len=3) :: sElem
    character(len=7) :: sCoeff
    integer          :: Coeff,I,LenRest,iEl
    !-------
    fOk= .true. 
    sRest= trim(sFormul)//"_"
    vStoik= 0
    Coeff=  0; Z=0; nDiv=1
    
    do
      
      if (sRest(1:1)=='_') exit      !end of formula
      
      LenRest=len_trim(sRest)
      I=SCAN(sRest,'(')              !find 1st occurrence of '('
      
      if(I==0) exit                  !no '(' found --> end of string
      
      if(I>1) then
        sElem=sRest(1:I-1)           !ElementName
        sRest=sRest(I+1:LenRest)     !rest of string
        I=SCAN(sRest,')')
        if(I>1) then
          sCoeff=sRest(1:I-1)        !Element Coeff
          sRest=sRest(I+1:LenRest)   !rest of string
        end if
      end if
      
      !call Str_Append(sElem,3)
      call WrdToInt(sCoeff,Coeff)
      
      select case(trim(sElem))
      case(""); exit                 !end of formula
      case("DIV") ; nDiv=Coeff ; cycle
      case("/")   ; nDiv=Coeff ; cycle
      case("E")   ; Z= Coeff   ; cycle
      case("+")   ; Z= Coeff   ; cycle
      case("-")   ; Z=-Coeff   ; cycle
      end select
      
      !if(DebFormula) write(fTrc,'(A)') trim(sElem)
      iEl=Name_Index(sElem,vNameEle)
      if(iEl==0) then
        fOk=.false.
        return    !element not found in input element list
      else
        vStoik(iEl)=Coeff
        !if(iFirst==0) iFirst=iEl !iFirst==0 checks whether iFirst already allocated
        !if(DebFormula) write(fTrc,'(A13,A3,A1,I3,A1,I3)') "El iEl Coeff ",sElem,T_,iEl,T_,Coeff  !,T_,Z
      end if
    
    end do
    
    return
  end subroutine Formula_Arxim_Read

end module M_Formula_Arxim
