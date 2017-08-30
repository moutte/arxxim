module M_Dtb_Read_Tools
!.tools for database reading routines
  use M_Kinds
  !
  implicit none
  !
  private
  !
  public:: DtbRead_Build_vElement
  public:: DtbRead_Build_ExtendedFormula
  !
contains

subroutine DtbRead_Build_vElement(vEle,vElement)
  use M_IOTools,only: cUpper,cLower
  use M_T_Element,only: T_Element
  !
  type(T_Element), intent(in) :: vEle(:)
  character(len=2),intent(out):: vElement(:)
  !
  character(len=2):: ss
  integer:: i
  !
  vElement= "  "
  do i= 1, size(vEle)
    ss="  "
    ss(1:1)= cUpper(vEle(i)%NamEl(1:1))
    if(vEle(i)%NamEl(2:2)/="_") ss(2:2)= cLower(vEle(i)%NamEl(2:2))
    vElement(i)= trim(ss)
  end do
  !
  vElement(size(vEle)+1)= '+'
  vElement(size(vEle)+2)= '-'
  vElement(size(vEle)+3)= '/'
  !
end subroutine DtbRead_Build_vElement

subroutine DtbRead_Build_ExtendedFormula(ff,vElement,W,isOk)
!--
!-- transform compact formula to extended formula
!--
  use M_Trace,  only: t_
  use M_Formula_Parser
  !
  integer,intent(in):: ff
  character(len=2),intent(in):: vElement(:)
  character(len=*),intent(inout):: W
  logical,intent(out):: isOk 
  !
  !! character(len=10):: FilCode
  !
  real(dp) :: b(size(vElement))
  real(dp) :: divreal
  character(len=71) :: s
  !
  if(ff>0) write(ff,'(A72,1X)',advance="NO") trim(W)
  !
  call Formula_Vector_Read ( W, vElement, s, b, isok, divreal )
  !
  if(.not. isok) then
    if(ff>0) write(ff,'(A)') " -> !!!CHECK!!!"
    return
  end if
  !
  !if(ff>0 .and. nint(divreal)/=1) &
  !& write(ff,'(A,I3,A)',advance="NO") "divreal=",nint(divreal)," = "
  !
  call Formula_Arxim_Build ( vElement, nint(b*divreal), W )
  !
  if(nint(divreal)/=1) then
    s=''
    write(s,'(I6)') nint(divreal)
    s= adjustl(s)
    W= trim(W)//'/('//trim(s)//')'
  end if
  !
  if(ff>0) write(ff,'(A)') trim(W)
end subroutine DtbRead_Build_ExtendedFormula

end module M_Dtb_Read_Tools
