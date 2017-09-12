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
  public:: FieldList_Read
  public:: FieldList_Check
  !
contains

subroutine FieldList_Read( &
& L,         &
& vStrField, &
& vifield)
!-- check that the essential fields are present in the field list:
!-- NAME, ECFORM/SCFORM, PARAMETERS 
  use M_Trace,   only: iDebug, Stop_
  use M_IoTools, only: LinToWrd
  !
  character(len=*),intent(inout):: L
  character(len=12),intent(in):: vStrField(:)
  integer,intent(out):: vifield(:)
  !
  character(len=80):: W
  logical:: EoL
  integer:: I,J
  !
  vifield(:)= 0
  I=0
  do
    call LinToWrd(L,W,EoL)
    I=I+1
    do J=1,size(vStrField)
      if( trim(W)==trim(vStrField(J)) ) then
        vifield(J)= I
        exit
      end if
    end do
    if(EoL) exit
  end do
  !
  if(iDebug==4) then
    do I=1,size(vStrField)
      print *,vifield(I),trim(vStrField(I))
    end do
  end if
  !
  if(vifield(10)==0) & ! for "PARAMETERS"
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(10)))
  !
  if(vifield(3)==0) & ! for "NAME"
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(3)))
  !
  if(vifield(4)==0 .and. vifield(5)==0) & ! for ECFORM/SCFORM
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(4))//"_"//trim(vStrField(5)))
  !
end subroutine FieldList_Read

subroutine FieldList_Check(i,vStrField,vifield)
!-- check that the essential fields are present in the field list:
!-- NAME, ECFORM/SCFORM, PARAMETERS 
  use M_Trace,   only: Stop_
  !
  integer,intent(in):: i
  character(len=12),intent(in):: vStrField(:)
  integer,intent(in):: vifield(:)
  !
  if(vifield(i)==0) &
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(i)))
  !
end subroutine FieldList_Check

subroutine DtbRead_Build_vElement(vEle,vElement)
!--
!-- transform uppercase element name to case-sensitive element name,
!-- e.g. MG to Mg, CO to Co, etc.
!-- necessary for scanning SCFORM formulas for stoichiometry retrieve
!--
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
