!! *******************************************************************
!! * File name:     M_echomat.f90
!! *
!! * Purpose:       - affichage controle de matrices, de vecteurs
!! *
!! * Author:        Anthony Michel (anthony.michel@ifp.fr)
!! *
!! * Created:       2005
!! *
!! * Modification:  2005
!! *
!! *
!! ********************************************************************

module M_echomat
  use m_kinds
  
  implicit none
  private

  character(len=8) :: Sfd='A' ! format affichage reel
  character(len=8) :: Rfd='EN20.10' ! format affichage reel
  !!character(len=8) :: Rfd='F24.10' ! format affichage reel
  character(len=8) :: Ifd='I6'    ! format affichage entiers

  public :: echomat_real
  public :: echomat_int
  public :: echovec_real
  public :: echovec_int
  public :: echovec_string

  public :: set_Rfd
  public :: set_Sfd
  public :: set_Ifd

contains

   subroutine set_Sfd(xSfd)

    implicit none
    character(len=*) :: xSfd

    Sfd = xSfd

  end subroutine set_Sfd

  !---------

  subroutine set_Rfd(xRfd)

    implicit none
    character(len=*) :: xRfd

    Rfd = xRfd

  end subroutine set_Rfd

  !---------

  subroutine set_Ifd(xIfd)

    implicit none
    character(len=*) :: xIfd

    Ifd = xIfd

  end subroutine set_Ifd

  !---------

  subroutine echomat_real(A, n, m, legend, trtag, unit)
    !------------------------------------------------------------
    !
    ! purpose : affichage de la matrice A
    !           si trtag='T', on transpose a l'affichage
    !
    !-----------------------------------------------------------

    implicit none

    real(dp), dimension(:,:) :: A
    integer :: i, j
    character(len=30) :: myform
    integer :: n, m
    character(len=*), optional :: legend
    character(len=*), optional :: trtag

    integer, optional :: unit
    integer :: iunit 
    !-------------
    iunit = 6
    if (present(unit)) iunit=unit

    if (present(legend)) then
       write(iunit,*) '----------- MATRIX : ', legend,'  ----------------'
    end if

    if (present(trtag).and.(trtag=='T')) then
       write(myform,'(1A1,1I3,2A)') '(',n,trim(Rfd),')'

       do i=1, m
          write(iunit,myform) (A(j,i), j=1, n)
       end do
    else
       write(myform,'(1A1,1I3,2A)') '(',m,trim(Rfd),')'

       do i=1, n
          write(iunit,myform) (A(i,j), j=1, m)
       end do
    end if

  end subroutine echomat_real

  !---------

  subroutine echomat_int(A, n, m, legend, trtag, unit)
    !------------------------------------------------------------
    !
    ! purpose : affichage de la matrice A
    !           si trtag='T', on transpose a l'affichage
    !
    !-----------------------------------------------------------

    implicit none

    integer, dimension(:,:) :: A
    integer :: i, j
    character(len=30) :: myform
    integer :: n, m
    character(len=*), optional :: legend
    character(len=*), optional :: trtag

    integer, optional :: unit
    integer :: iunit = 6
    !-------------
    if (present(unit)) iunit=unit

    if (present(legend)) then
       write(iunit,*) '----------- MATRIX : ', legend,'  ----------------'
    end if

    if (present(trtag).and.(trtag=='T')) then
       write(myform,'(1A1,1I3,2A)') '(',n,trim(Ifd),')'
       do i=1, m
          write(iunit,myform) (A(j,i), j=1, n)
       end do
    else
       write(myform,'(1A1,1I3,2A)') '(',m,trim(Ifd),')'
       do i=1, n
          write(iunit,myform) (A(i,j), j=1, m)
       end do
    end if

  end subroutine echomat_int

  !---------

  subroutine echovec_int(V, n, legend, trtag, unit)
    !------------------------------------------------------------
    !
    ! purpose : affichage de la matrice A
    !           si trtag='T', on transpose a l'affichage
    !
    !-----------------------------------------------------------

    implicit none

    integer, dimension(:) :: V
    integer :: i, j
    character(len=30) :: myform
    integer :: n !, m
    character(len=*), optional :: legend
    character(len=*), optional :: trtag

    integer, optional :: unit
    integer :: iunit = 6
    !-------------
    if (present(unit)) iunit=unit

    if (present(legend)) then
       write(iunit,*) '----------- VECTOR : ', legend,'  ----------------'
    end if

    if (present(trtag).and.(trtag=='T')) then
       write(myform,'(1A1,1I3,2A)') '(',n,trim(Ifd),')'
       write(iunit,myform) (V(j), j=1, n)
    else
       write(myform,'(1A1,1I3,2A)') '(',1,trim(Ifd),')'
       do i=1, n
          write(iunit,myform) V(i)
       end do
    end if

  end subroutine echovec_int

!-----


  subroutine echovec_string(V, n, legend, trtag, unit)
    !------------------------------------------------------------
    !
    ! purpose : affichage de la matrice A
    !           si trtag='T', on transpose a l'affichage
    !
    !-----------------------------------------------------------

    implicit none

    character(len=*), dimension(:) :: V
    integer :: i, j
    character(len=30) :: myform
    integer :: n
    character(len=*), optional :: legend
    character(len=1), optional :: trtag

    integer, optional :: unit
    integer :: iunit = 6
    logical :: LT
    !-------------
    LT=.false.

    if (present(unit)) iunit=unit

    if (present(legend)) then
       write(iunit,*) '----------- VECTOR : ', legend,'  ----------------'
    end if

    if (present(trtag)) then
       if (trtag=='T')  LT=.true.
    end if

    if (LT) then
       write(myform,'(1A1,1I3,2A)') '(',n,trim(Sfd),')'
       write(iunit,myform) (trim(V(j))//'  ', j=1, n)
    else
       write(myform,'(1A1,1I3,2A)') '(',1,trim(Sfd),')'
       do i=1, n
          write(iunit,myform) trim(V(i))
       end do
    end if

  end subroutine echovec_string


  !-----


  subroutine echovec_real(V, n, legend, trtag, unit)
    !------------------------------------------------------------
    !
    ! purpose : affichage de la matrice A
    !           si trtag='T', on transpose a l'affichage
    !
    !-----------------------------------------------------------

    implicit none

    real(dp), dimension(:), intent(in) :: V
    integer :: i, j
    character(len=60) :: myform
    integer, intent(in) :: n
    character(len=*), optional :: legend
    character(len=1), optional :: trtag

    integer, optional :: unit
    integer :: iunit
    logical :: LT
    !-------------
    LT = .false.
    iunit = 6
    
    if (present(unit)) iunit=unit
        
    if (present(legend)) then
       write(iunit,*) '----------- VECTOR : ', legend,'  ----------------'
    end if

    if (present(trtag)) then
       if (trtag=='T')  LT=.true.
    end if

    if (LT) then
       write(myform,'(1A1,1I3,2A)') '(',n,trim(Rfd),')'
       write(iunit,myform) (V(j), j=1, n)
    else
       write(myform,'(1A1,1I3,2A)') '(',1,trim(Rfd),')'
       do i=1, n
          write(iunit,myform) V(i)
       end do
    end if

  end subroutine echovec_real


end module M_echomat
