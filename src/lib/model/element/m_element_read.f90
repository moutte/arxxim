module M_Element_Read
  use M_Kinds
  use M_Trace,    only: iDebug,T_,fTrc,Stop_,Warning_
  use M_T_Element,only: T_Element
  implicit none

  private
  
  public:: Elements_BuildLnk
  public:: Elements_LnkToVec
  public:: LnkEle_Build
  public:: T_LnkEle
  public:: Element_Read_Redox
  public:: Element_Read_Entropy
  !
  type:: T_LnkEle
    type(T_Element)         ::Value
    type(T_LnkEle), pointer ::Next
  end type T_LnkEle
  !
contains

subroutine Element_Read_Entropy(vEle)
!--
!-- read entropy of elements at standard state
!--
!-- entropy data for elements
!-- are needed for conversion from Berman-Brown
!-- to Benson-Helgeson convention
!--
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  use M_Files,    only: NamFInn
  use M_T_Element,only: T_Element,Element_Index
  !
  type(T_Element),intent(inout):: vEle(:)
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios,i
  !
  if(idebug>1) write(fTrc,'(/,A,/)') "< Elements_Read_Entropy"
  !
  call GetUnit(F)
  open(F,file=trim(NamFInn))
  !
  DoFile: do 
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToend(L,W,EoL)
    !
    select case(W)
    case("ELEMENTS.ENTROPY") !!!!!!canevas for reading one "block"
      DoBlock: do
        read(F,'(A)',iostat=ios) L
        if(ios/=0)                                         exit DoFile
        call LinToWrd(L,W,EoL)
        call AppendToend(L,W,EoL)
        if(W(1:1)=='!')                                    cycle DoBlock
        select case(W)
        case("END","ENDELEMENTS.ENTROPY") ;                exit DoFile
        end select
        call Str_Append(W,3)
        i= Element_Index(W(1:3),vEle)
        if(i>0) then
          call LinToWrd(L,W,EoL)
          call WrdToReal(W,vEle(i)%S0)
        else
          call Warning_("Element "//trim(W)//" not in base")
        end if
      end do DoBlock
    end select
    !
  end do DoFile!  close(F)
  !----------------------------------------------------------------debug
  if(idebug>1) then
    write(fTrc,'(A,/)') "entropy of element at standard conditions:"
    do i=1,size(vEle)
      write(fTrc,"(A3,A1,G12.3)") &
      & vEle(I)%NamEl, T_, vEle(I)%S0
    end do
  end if
  !---------------------------------------------------------------/debug
  if(idebug>1) write(fTrc,'(/,A,/)') "</ Elements_Read_Entropy"
  !
end subroutine Element_Read_Entropy
  
subroutine Element_Read_Redox(vEle)
!--
!-- read redox state of (some) elements
!--
!-- format of block:
!--   ELEMENTS.REDOX
!--     FE 2
!--     ../..
!--   end
!--
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  use M_Files,    only: NamFInn
  use M_T_Element,only: T_Element,Element_Index
  !
  type(T_Element),intent(inout):: vEle(:)
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios,i
  !
  if(idebug>1) write(fTrc,'(/,A,/)') "< Elements_Read_Redox"
  !
  call GetUnit(F)
  open(F,file=trim(NamFInn))
  !
  DoFile: do 
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToend(L,W,EoL)
    select case(W)
    case("ELEMENTS.REDOX") !!!!!!canevas for reading one "block"
      DoBlock: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        call AppendToend(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlock !skip comment lines
        select case(W)
        case("END","ENDELEMENTS.REDOX"); exit DoFile
        end select
        !
        call Str_Append(W,3)
        i= Element_Index(W(1:3),vEle)
        if(i>0) then
          call LinToWrd(L,W,EoL)
          call WrdToInt(W,vEle(i)%Z)
        else
          call Warning_("Element "//trim(W)//" not in base")
        end if
        !
      end do DoBlock
    end select
  end do DoFile
  close(F)
  !
  if(idebug>1) then
    write(fTrc,'(A,/)') "oxydation states, default values:"
    do i=1,size(vEle)
      write(fTrc,"(A3,A1,I3)") &
      & vEle(I)%NamEl, T_, vEle(I)%Z
    end do
  end if
  !
  if(idebug>1) write(fTrc,'(/,A,/)') "</ Elements_Read_Redox"
end subroutine Element_Read_Redox
  
subroutine LnkEle_Build(B,E,L,P)
  logical        ::b
  type(T_Element)::E
  type(T_LnkEle),pointer::L,P
  !
  if(B) then
    nullify(L)
    allocate(L)
    nullify(L%next)
    L%Value=E
    P=>L
  else
    allocate(P%next)
    nullify(P%next%next)
    P%next%Value=E
    P=>P%next
  end if
end subroutine LnkEle_Build

subroutine Elements_LnkToVec(LnkEle,vEle)
  use M_T_Element,only: T_Element
  !
  type(T_LnkEle), pointer    :: LnkEle
  type(T_Element),intent(out):: vEle(:)
  !
  type(T_LnkEle),pointer:: pCur, pPrev
  integer:: I
  !
  I=0
  pCur=>LnkEle
  do while (associated(pCur))
    I= I+1
    vEle(I)=pCur%Value
    pPrev=>pCur
    pCur=> pCur%next
    deallocate(pPrev)
  end do
  !
end subroutine Elements_LnkToVec

subroutine Elements_BuildLnk(LnkEle,nEle)
  use M_IOTools !,only: LinToWrd
  use M_Files,    only: NamFEle
  use M_T_Element,only: T_Element
  !
  type(T_LnkEle),pointer    :: LnkEle
  integer,       intent(out):: nEle
  !
  logical:: L_Aqueous=.true.
  !
  type(T_Element)   :: Ele
  character(len=512):: L,W,sListElem
  logical           :: EoL
  integer           :: f,ios
  type(T_LnkEle),pointer::pCur
  !
  if(idebug>1) write(fTrc,'(A)') "< Elements_BuildLnk"
  !
  ! nEle-> total nr in database, not the nCp of the run
  !
  call GetUnit(f)
  open(f,file=trim(NamFele))
  !
  !--------------------------------------------------- build linked list
  DoFile: do
    read(F,'(A)',iostat=ios) L  ; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=="!") cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    select case(W)
    !
    case("ENDINPUT")  ; exit DoFile
    !
    case("ELEMENTS")
      nEle=0
      sListElem=""
      !
      if(L_Aqueous) then
        Ele%NamEl="O__"  ;  Ele%WeitKg=0.015999D0  ;  Ele%Z=-2  ;  Ele%S0=102.57D0
        nEle=nEle+1
        call LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
        sListElem=trim(sListElem)//trim(Ele%NamEl)
        !
        Ele%NamEl="H__"  ;  Ele%WeitKg=0.001008D0  ;  Ele%Z=1   ;  Ele%S0=65.34D0
        nEle=nEle+1
        call LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
        sListElem=trim(sListElem)//trim(Ele%NamEl)
        !
        Ele%NamEl="OX_"  ;  Ele%WeitKg=Zero        ;  Ele%Z=1   ;  Ele%S0=Zero
        nEle=nEle+1
        call LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
        sListElem=trim(sListElem)//trim(Ele%NamEl)
        !
        !-> sListElem="O__H__OX_"
        if(iDebug==5) print '(A)',trim(sListElem)
      end if
      !
      LoopReadElem: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=="!") cycle LoopReadElem !skip comment lines
        call AppendToEnd(L,W,EoL)
        select case(W)
          case("ENDINPUT"); exit DoFile
          case("END","ENDELEMENTS"); exit LoopReadElem
        end select
        
        call Str_Append(W,3)
        if(index(sListElem,trim(W))<1) then
          !
          call Str_Append(W,3)
          Ele%NamEl=trim(W)
          if(idebug>1) write(fTrc,"(A3)") Ele%NamEl
          !
          sListElem=trim(sListElem)//trim(W)
          if(iDebug==5) print '(A)',trim(sListElem)
          !
          call LinToWrd(L,W,EoL); call WrdToReal(W,Ele%WeitKg)
          !
          call LinToWrd(L,W,EoL); call WrdToInt (W,Ele%Z)
          if(Ele%Z==0) then  ; Ele%Redox="VAR"
          else               ; Ele%Redox="FIX"
          end if
          call LinToWrd(L,W,EoL)
          !
          call LinToWrd(L,W,EoL); call WrdToReal(W,Ele%S0)
          !
          nEle=nEle+1
          !!if(nEle <= nElMax) then
          call LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
          !!else
          !!  call Warning_("Too many elements >> element Skipped")
          !!end if
        end if
        !
      end do LoopReadElem
      
      if(idebug>1) write(fTrc,"(A,I3)") "nEle=",nEle
      
    !endcase("ELEMENT")
    end select
  end do DoFile
  close(f)
  !--------------------------------------------------/ build linked list
  !
  if(nEle==0) call Stop_("Found NO Elements ... (missing path ??)")
  !
  if(idebug>1) write(fTrc,'(A)') "</ Elements_BuildLnk"
end subroutine Elements_BuildLnk

end module M_Element_Read
