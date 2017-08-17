module M_T_Element
  use M_Kinds
  !
  implicit none
  !
  private
  !-- note on structured data types used in the project --
  !
  !-- at the moment, all data types implemented are of fixed size, 
  !-- they do not contain fields of allocatable or pointer type,
  !
  !-- links between objects are realized using indexes,
  !-- e.g. a component c comprises an integer index C%iEle 
  !-- that points to the element vEle(C%iEle) of the element list vEle
  !
  !-- character strings are all of fixed size
  !-- (generally - n*8 -1, for better output on text editor using tabs of length 8 ...)
  !
  public:: T_Element
  public:: Element_Index
  public:: Element_Zero
  public:: Formula_Read
  public:: Formula_Build
  !
  type:: T_Element
  !-- container for storing a CHEMICAL ELEMENT
    character(len=3):: NamEl
    character(len=3):: Redox !either "FIX", or "VAR" (for Cl, Fe, Cr, As, ...)
    !=> whether the default oxydation state of the element is fixed or variable !!
    real(dp)     :: &
    & WeitKg, & ! atomic weight (kg)
    & S0        ! entropy (J) at 25°C/1atm
    !             (used for conversions between thermodyn. conventions)
    integer      :: Z 
    ! Z= nominal valency of element, used for computing species charge, redox state, etc.
  end type T_Element
  type(T_Element):: Element_Zero= T_Element("ZZZ","VAR",One,Zero,0)
  !
contains  
 
!! type(T_Element) function Element_Zero
!!   Element_Zero%Name=   "Z"
!!   Element_Zero%Redox=  "VAR"
!!   Element_Zero%WeitKg= One
!!   Element_Zero%S0=     Zero
!!   Element_Zero%Z=      0
!! end function Element_Zero

integer function Element_Index(Str,V)
!--
!-- -> position of T_element with %Name--Str in V
!--
  character(len=*),intent(in):: Str
  type(T_Element), intent(in):: V(:)
  integer     ::I
  !
  Element_Index=0
  if(size(V)==0) return
  !
  I=0
  do
    I=I+1 
    if(trim(Str)==trim(V(I)%NamEl)) then
      Element_Index=I
      exit
    end if
    if(I==size(V)) exit
  end do !if Str not found -> I=0
  !
  return
end function Element_Index

!---

subroutine Formula_Read( &
& sFormul,vEle, &
& Z,nDiv,fOk,vStoik)

  use M_Formula_Parser
  
  implicit none
  !
  character(len=*),intent(in) :: sFormul !the stoichio'formula of a substance
  type(T_Element), intent(in) :: vEle(:) !the current element list
  !
  integer,intent(out):: Z,nDiv ! the charge and the "divisor"
  logical,intent(out):: fOk    ! the formula is consistent (with elemental basis)
  integer,intent(out):: vStoik(size(vEle)) !the "formula vector"
  
  call Formula_Arxim_Read_Standard (sFormul,vEle(:)%NamEl,Z,nDiv,fOk,vStoik)
  
  return
end subroutine Formula_Read

!---

subroutine Formula_Build(vEle,vStoik,Zsp,nDiv,S) !,Discret_,nDiv1,nDiv2,I_,J_,S)
  use M_Formula_Parser
  implicit none
  !build a formula with a "static" element order, -> to make species sorting easier ...
  type(T_Element), intent(in) :: vEle(:)
  integer,         intent(in) :: vStoik(:)
  integer,         intent(in) :: Zsp,nDiv
  character(len=*),intent(out):: S

  call Formula_Arxim_Build_Standard(vEle(:)%NamEl,vStoik,Zsp,nDiv,S)
  
end subroutine Formula_Build

end module M_T_Element

