!! *******************************************************************
!! * File name:     M_iolog.f90
!! *
!! * Purpose:       tracer les ouvertures et fermetures de fichiers
!!
!! * Author:        Anthony Michel (anthony.michel@ifp.fr)
!! *
!! * Created:       2006
!! *
!! * Modification:  2006
!! *
!! *
!! ********************************************************************

module M_iolog

  implicit none
  private
  
  public :: iolog_open   ! tri ordre lexicographique
  public :: iolog_close  ! tri ordre lexicographique
  public :: init_iolog
  public :: end_iolog
  integer :: iolog_unit = 3000
  integer :: nbfile = 0
  logical, public :: active_iolog = .false.

contains

  subroutine init_iolog(filename)

    implicit none
    character(len=*), intent(in) :: filename
    !--------
    nbfile =0
    if (active_iolog) then
       open(unit=iolog_unit, file=filename, status='UNKNOWN')
       write(iolog_unit, '(A,I3,A, A, A, A)') &
            '[',nbfile ,']', '[ ', 'INIT',' ]'
       call iolog_open(unit=iolog_unit, file=filename, status='UNKNOWN')
    end if

  end subroutine init_iolog

  !---

  subroutine end_iolog

    implicit none

    if (active_iolog) then
       call iolog_close(unit=iolog_unit)
       write(iolog_unit, '(A,I3,A, A, A, A)') &
            '[',nbfile ,']', '[ ', 'end',' ]'
       close(unit=iolog_unit)
    end if
  end subroutine end_iolog
  
  !---

  subroutine iolog_open(unit, file, status)
    
    implicit none
    character(len=*), intent(in) :: status
    character(len=*), intent(in) :: file
    integer, intent(in) :: unit
    character(len=10) :: lstatus
    !--------
   
    nbfile=nbfile+1
    if (active_iolog) then
       if (status=='') then
          lstatus = 'UNKNOWN'
       else
          lstatus = status
       end if
           
       write(iolog_unit, '(A, I3, A, A, A, A, I10, A, A30,A, A8, A)') &
            '[',nbfile ,']','[ ', 'open',' | ', unit , ' | ', trim(adjustl(file)), ' | ', lstatus, ' ]'

    end if
  end subroutine iolog_open

  !---

  subroutine iolog_close(unit)

    implicit none
    integer, intent(in) :: unit
    !--------
    nbfile = nbfile - 1
    if (active_iolog) then
       write(iolog_unit, '(A, I3, A, A, A,A,I10,A)') &
            '[',nbfile ,']','[ ', 'close',' | ', unit , ' ]'
    end if
  end subroutine iolog_close
  

end module M_iolog
