 module M_Formula_Arxim_Standard
   
  !===============================================================
  ! Purpose :  Arxim Formulas Utils Standard
  !---------------------------------------------------------------
  ! Arxim Formula Syntax 
  !---------------------------------------------------------------
  ! CaCO3        = Ca(3)C(1)O(3) 
  ! HCO3-        = H(1)C(1)O(3)-(1)
  ! Ca+2         = Ca(1)+(2)
  ! CO3(MgCa)0.5 = C(2)O(6)MG(1)CA(1)/(2) 
  !             or C(2)O(6)MG(1)CA(1)DIV(2)
  !===============================================================
  implicit none
  private

  public :: Formula_Arxim_Build_Standard
  public :: Formula_Arxim_Read_Standard

  logical, parameter :: DebFormula = .false.

  !---

contains
  
  
  subroutine Formula_Arxim_Read_Standard(sFormul,vNameEle,Z,nDiv,fOk,vStoik) 
  !=============================================================================
  !-- read a chemical formula sFormul
  !-- to a stoikiometric vector (vStoik)/Div and the charge Z
  !-- according to the element list vEle
  !--
  !---> used to read the elemental decompositon of a substance from its formula
  !--
  !-- CAVEAT
  !--   vStoik is an array of integers; nDiv the formula divisor
  !--   thus, the "real" stoichio coeff for element vEle(i) is vstoik(i)/Div !!!
  !--
  !=============================================================================

    use M_Formula_Utils
    use M_IOTools,only: Str_Append, WrdToInt
    implicit none
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
    
    if(trim(sRest)=="_") then
      fOk= .false.
      return
    end if
    
    vStoik= 0
    Coeff=  0
    Z=      0
    nDiv=   1
    
    do
      if (sRest(1:1)=='_') exit               !end of formula
      
      LenRest= len_trim(sRest)
      I= SCAN(sRest,'(')                     !find 1st occurrence of '('
      
      if(I==0) exit                          !no '(' found --> end of string
      
      if(I>1) then
      
        sElem=sRest(1:I-1)                   !ElementName
        sRest=sRest(I+1:LenRest)             !rest of string
        I=SCAN(sRest,')')
        if(I>1) then
          sCoeff=sRest(1:I-1)                !Element Coeff
          sRest=sRest(I+1:LenRest)           !rest of string
        end if
      
      end if
       
      call Str_Append(sElem,3)
      call WrdToInt(sCoeff,Coeff)
       
      select case(trim(sElem))
        case("___")  ; exit                    !end of formula
        case("DIV")  ; nDiv=Coeff ;   cycle
        case("/__")  ; nDiv=Coeff ;   cycle
        case("E__")  ; Z= Coeff   ;   cycle
        case("+__")  ; Z= Coeff   ;   cycle
        case("-__")  ; Z=-Coeff   ;   cycle
      end select
       
      iEl=Name_Index(sElem,vNameEle)
       
      if(iEl==0) then
      
        fOk=.false.  !element not found in input element list
        return
       
      else
      
        vStoik(iEl)=Coeff
        ! if(iFirst==0) iFirst=iEl !iFirst==0 checks whether iFirst already allocated
        ! if(DebFormula) write(fTrc,'(A13,A3,A1,I3,A1,I3)') &
        ! & "El iEl Coeff ",sElem,T_,iEl,T_,Coeff  !,T_,Z
        
      end if
       
    end do
    
    return
  end subroutine Formula_Arxim_Read_Standard

  !---
  
  subroutine Formula_Arxim_Build_Standard(vNameEle,vStoik,Zsp,nDiv,S) !,Discret_,nDiv1,nDiv2,I_,J_,S)
    !--===============================================
    ! Purpose : build a formula with a "static" element order,
    !           -> to make species sorting easier ...
    !--===============================================
    implicit none
    character(len=*), intent(in) :: vNameEle(:)
    integer,          intent(in) :: vStoik(:)
    integer,          intent(in) :: Zsp,nDiv
    character(len=*), intent(out):: S
    !
    character(len=4):: Str
    integer:: iEl
    S=""
    do iEl=size(vNameEle),1,-1 !-> will have H and O as last elements 
      if(vStoik(iEl)>0) then
        write(Str,'(I4)') vStoik(iEl) !; if(vStoik(iEl)<10) Str(1:1)='0'
        S= trim(S)//vNameEle(iEl)(1:2)//"("//trim(adjustl(Str))//")"
      end if
    end do
    if(Zsp>0) then
      write(Str,'(I4)') Zsp  ; S= trim(S)//"+("//trim(adjustl(Str))//")"
    end if
    if(Zsp<0) then
      write(Str,'(I4)') -Zsp  ; S= trim(S)//"-("//trim(adjustl(Str))//")"
    end if
    if(nDiv>1) then
      write(Str,'(I4)') nDiv  ; S= trim(S)//"/("//trim(adjustl(Str))//")"
    end if
  end subroutine Formula_Arxim_Build_Standard

end module M_Formula_Arxim_Standard
