module M_Formula_Utils
  implicit none
  
  private
  
  public :: Name_Index

  contains
  
  !---
  
  function Name_Index(Str,V) result(index)
    !=====================================================
    ! Position of Name with Name==Str in VName
    !=====================================================
    use m_trace, only : iDebug
    implicit none
    character(len=*),              intent(in)::Str
    character(len=*), dimension(:),intent(in)::V
    !---
    integer :: I
    integer :: index
    !---
    index=0
    if(size(V)==0) return
    I=0
    do
      I=I+1 
      if(iDebug==5) write(*,*) "Testing" , trim(Str), trim(V(I))
      if(trim(Str)==trim(V(I))) then
        index=I
        exit
      end if
      if(I==size(V)) exit
    end do !if Str not found -> I=0
    return
  end function Name_Index
 
end module M_Formula_Utils
