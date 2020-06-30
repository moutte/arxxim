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
    character(len=7):: SysName= "ZZZ" ! MASTER, INJECT, BOX, ...
    character(len=7):: SysType= "ZZZ" ! AQUEOUS, ...
    !
    real(dp):: SysTdgK= 25.D0
    real(dp):: SysPbar= 1.D0
    !
    integer :: SysNCp= 0
    type(T_Component),allocatable:: vCpn(:)
    !
  end type T_System
  
contains

subroutine System_New(N,This)
  integer,intent(in):: N
  type(T_System),intent(out):: This
  !
  This%SysName= "ZZZ"
  This%SysType= "ZZZ"
  !
  This%SysTdgK= 25.D0
  This%SysPbar= 1.D0
  
  This%SysNCp=  N
  allocate(This%vCpn(N))
  !
end subroutine System_New

subroutine System_Free(This)
  type(T_System),intent(inout):: This
  !
  This%SysName="ZZZ"
  This%SysType="ZZZ"
  This%SysNCp= 0
  deallocate(This%vCpn)
  !
end subroutine System_Free

end module M_T_System

