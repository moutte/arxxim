module M_Box_System_Vars

  !------------------------------------------------------------
  ! Structure Parameters for Module Box
  !------------------------------------------------------------

  use M_Kinds
  use M_Trace
  
  use M_Box_Param_Vars, only: nAq, nMk, nCp, nPr, nAs
  use M_T_Component

  implicit none
  public

  !// System Size Parameters
  public :: nAq
  public :: nMk
  !--
  public :: nCp
  public :: nPr
  public :: nAs

  !// Components Vector
  type(T_Component), allocatable :: vCpnBox(:)

  !// Formula Matrix
  real(dp), allocatable:: tAlfAq(:,:)   ! 1:ncp, 1:nAq
  real(dp), allocatable:: tAlfMk(:,:)   ! 1:nCp, 1:nMk

  real(dp), allocatable:: tAlfPr(:,:)   ! 1:ncp, 1:nPr
  real(dp), allocatable:: tAlfAs(:,:)   ! 1:ncp, 1:nAs

  !// Stoikiometry
  real(dp), allocatable:: tNuAs(:,:)    ! 1:nAs, 1:nPr
  real(dp), allocatable:: tNuMk(:,:)    ! 1:nMk, 1:nPr

  !// Pointers to aqueous species
  integer, allocatable :: vOrdAq(:)
  integer :: iAqH2O
 
  !// Pointers to kinetic mineral species
  integer, allocatable :: vOrdMk(:)

contains

  subroutine Box_System_Vars_New

    implicit none

    call Info_("Box_System_Vars_New")

    !// Components Vector
    allocate(vCpnBox(nCp))

    !// Formula Matrix
    allocate(tAlfAq(nCp, nAq))
    allocate(tAlfMk(nCp, nMk))

    allocate(tAlfPr(nCp, nPr))
    allocate(tAlfAs(nCp, nAs))

    !// Stoikiometry
    allocate(tNuAs(nAs, nPr))
    allocate(tNuMk(nMk, nPr))    

    !// Pointers to aqueous species
    allocate(vOrdAq(nAq))

    !// Pointers to kinetic mineral species
    allocate(vOrdMk(nMk))

    !// Clean variables after allocation
    call Box_System_Vars_Zero

  end subroutine Box_System_Vars_New

  !---

  subroutine Box_System_Vars_Zero

    implicit none
    integer :: iCp

    call Info_("Box_System_Vars_Zero")

    !// Components Vector
    do iCp = 1, nCp
       call Component_Zero(vCpnBox(iCp))
    end do
    
    !// Formula Matrix
    tAlfAq(:,:)= Zero
    tAlfMk(:,:)= Zero

    tAlfAs(:,:)= Zero
    tAlfPr(:,:)= Zero

    !// Stoikiometry
    tNuAs(:,:)= Zero
    tNuMk(:,:)= Zero

    !// Pointers to aqueous species
    vOrdAq(:) = 0
    iAqH2O= 0
    
    !// Pointers to kinetic mineral species
    vOrdMk(:) = 0

  end subroutine Box_System_Vars_Zero

  !---

  subroutine Box_System_Vars_Delete

    implicit none

    call Info_("Box_System_Vars_Delete")

    !// Clean variables before deallocation
    call Box_System_Vars_Zero

    
    !// Components Vector
    deallocate(vCpnBox)

    !// Formula Matrix
    deallocate(tAlfAq)
    deallocate(tAlfMk)

    deallocate(tAlfPr)
    deallocate(tAlfAs)

    !// Stoikiometry
    deallocate(tNuAs)
    deallocate(tNuMk)

    !// Pointers to aqueous species
    deallocate(vOrdAq)

    !// Pointers to kinetic mineral species
    deallocate(vOrdMk)

  end subroutine Box_System_Vars_Delete

end module M_Box_System_Vars
