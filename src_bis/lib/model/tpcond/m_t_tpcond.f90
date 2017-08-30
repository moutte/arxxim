module M_T_TPCond
!--
!-- T,P conditions
!--

  use M_Kinds
  implicit none
  !
  private
  !
  public:: T_TPCond
  !
  !! public:: TPcond_IndexTP
  !! public:: TPcond_IndexT
  !
  type:: T_TPCond !container for T,P Properties
    character(15)::Name
    real(dp):: TdgC, Pbar
  end type T_TPCond
  
  type:: T_TPSeries !container for a series of T,P Properties
    character(15)::Name
    integer:: NPoints
    real(dp),allocatable:: vTdgC(:), vPbar(:)
  end type T_TPSeries
  
  real(dp),parameter:: Delta=1.0D0

contains

!!  if(allocated(vTPpath)) deallocate(vTPpath)

integer function TPcond_IndexTP(TdgK,Pbar,vTPCond)
  use M_Dtb_Const,only: T_CK
  real(dp),       intent(in):: TdgK,Pbar
  type(T_TPCond), intent(in):: vTPCond(:)
  integer :: I,J
  J=0
  do I=1,size(vTPCond)
    if(   ABS(TdgK -T_CK-vTPCond(I)%TdgC)<Delta &
    .and. ABS(Pbar      -vTPCond(I)%Pbar)<Delta ) J=I
  enddo
  TPcond_IndexTP= J
end function TPcond_IndexTP

integer function TPcond_IndexT(TdgK,vTPCond)
  use M_Dtb_Const,only: T_CK
  real(dp),      intent(in):: TdgK
  type(T_TPCond),intent(in):: vTPCond(:)
  integer :: I,J
  J=0
  do I=1,size(vTPCond)
    if(ABS(TdgK -T_CK-vTPCond(I)%TdgC)<Delta) J=I
  enddo
  TPcond_IndexT= J
end function TPcond_IndexT

end module M_T_Tpcond

