module M_Clean_Jac
  !-----------------------------------------------------------
  ! Clean Jacobian Matrix from Small Coupling terms
  ! Cleaning operation based on max absolute value term
  !-----------------------------------------------------------
  use M_Kinds
  use M_Trace
  implicit none
  private

  real(dp), parameter :: xeps = 1D-10
  
  public :: Clean_Jac_All
  public :: Clean_Jac_Index
  public :: Clean_Jac_Section
  
contains 

  subroutine Clean_Jac_All(tJac)
    !-------------------------------------------------
    ! Clean All Equations
    !-------------------------------------------------
    implicit none
    real(dp) :: tJac(:,:)
    integer :: nEq 
    !--
    !call Info_( "Clean_Jac_All" )
    nEq = size(tJac,1)
    call Clean_Jac_Section(tJac, 1, nEq)
   
  end subroutine Clean_Jac_All

  !---
  
  subroutine Clean_Jac_Section(tJac, nA, nB)
    !--------------------------------------------------------
    ! Clean Equations [nA:nB] 
    !--------------------------------------------------------
    implicit none
    real(dp) :: tJac(:,:)
    integer :: nA, nB
    real(dp) :: x, xmax
    integer :: iEq,j, N
    !
    !call Info_( "Clean_Jac_Section" )
    N = size(tJac,2)
    do iEq=nA,nB
       ! compute maximal species contribution for component iC
       xmax = maxval ( abs( tJac(iEq,1:N)) ) 
       do j=1,N
          ! clean small relative contributions of the line
          x = tJac(iEq,j)
          !write(*,*) "x=", x, xmax
          if (abs(x)<xeps*xmax) then
             !write(*,*) "clean value tJac"
             tJac(iEq,j) = Zero
          end if
       end do
    end do

  end subroutine Clean_Jac_Section

  !--
  
  subroutine Clean_Jac_Index(tJac, IdxEq)
    !--------------------------------------------------------
    ! Clean Indexed Equations
    !--------------------------------------------------------
    implicit none
    real(dp) :: tJac(:,:)
    integer  :: IdxEq(:)
    !---
    real(dp) :: x, xmax
    integer  :: iEq,j
    integer  :: N, nEq
    !call Info_( "Clean_Jac_Index" )
    N   = size(tJac,2)
    nEq = size(IdxEq)
    
    do iEq = IdxEq(1),IdxEq(nEq)
       ! compute maximal species contribution for component iC
       xmax = maxval ( abs( tJac(iEq,1:N)) ) 
       do j=1,N
          ! clean small relative contributions of the line
          x = tJac(iEq,j)
          if (abs(x)<xeps*xmax) tJac(iEq,j) = Zero
       end do
    end do

  end subroutine Clean_Jac_Index
  
end module M_Clean_Jac
