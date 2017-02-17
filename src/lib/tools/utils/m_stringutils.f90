module M_StringUtils
  
  implicit none
  public

contains

  integer function String_Index(vNames, Str) !,nSpc_)
    !.Position of names in the Names Tableau
    character(len=*),             intent(in):: Str
    character(len=*),dimension(:),intent(in):: vNames
    !
    integer     ::I, N
    String_Index=0
    N = size(vNames)
    if(N==0) return
    I=0
    do I=1,N
       if(trim(Str)==trim(vNames(I))) then
          String_Index=I
          exit
       end if
    end do
    return
  end function String_Index

end module M_StringUtils
