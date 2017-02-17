module M_T_System

  use M_Kinds
  use M_T_Component,only: T_Component
  
  implicit none
  
  private
  
  public:: &
  & T_System,   &
  & System_New, &
  & System_Free
  
  type:: T_System
    !private
    character(len=7):: SysName ! MASTER, INJECT, BOX, ...
    character(len=7):: SysType ! AQUEOUS, ...
    !
    real(dp):: SysTdgK, SysPbar
    !
    integer :: SysNCp
    type(T_Component),allocatable:: vCpn(:)
    !
  end type T_System
  
contains

subroutine System_New(N,This)
  integer,intent(in):: N
  type(T_System),intent(out):: This
  !
  This%SysName= "VOID"
  This%SysType= "VOID"
  
  This%SysTdgK= 25.D0
  This%SysPbar= 1.D0
  
  This%SysNCp=  N
  allocate(This%vCpn(N))
  !
end subroutine System_New

subroutine System_Free(This)
  type(T_System),intent(inout):: This
  !
  This%SysName="VOID"
  This%SysType="VOID"
  This%SysNCp= 0
  deallocate(This%vCpn)
  !
end subroutine System_Free

end module M_T_System

