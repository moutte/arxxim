MODULE M_Clean_Jac
  !-----------------------------------------------------------
  ! Clean Jacobian Matrix from Small Coupling terms
  ! Cleaning operation based on max absolute value term
  !-----------------------------------------------------------
  USE M_Kinds
  USE M_Trace
  IMPLICIT NONE
  PRIVATE

  REAL(kind=8), PARAMETER :: xeps = 1D-10
  
  PUBLIC :: Clean_Jac_All
  PUBLIC :: Clean_Jac_Index
  PUBLIC :: Clean_Jac_Section
  
CONTAINS 

  SUBROUTINE Clean_Jac_All(tJac)
    !-------------------------------------------------
    ! Clean All Equations
    !-------------------------------------------------
    IMPLICIT NONE
    REAL(dp) :: tJac(:,:)
    INTEGER :: nEq 
    !--
    !CALL Info_( "Clean_Jac_All" )
    nEq = SIZE(tJac,1)
    CALL Clean_Jac_Section(tJac, 1, nEq)
   
  END SUBROUTINE Clean_Jac_All

  !---
  
  SUBROUTINE Clean_Jac_Section(tJac, nA, nB)
    !--------------------------------------------------------
    ! Clean Equations [nA:nB] 
    !--------------------------------------------------------
    IMPLICIT NONE
    REAL(dp) :: tJac(:,:)
    INTEGER :: nA, nB
    REAL(kind=8) :: x, xmax
    INTEGER :: iEq,j, N
    !
    !CALL Info_( "Clean_Jac_Section" )
    N = SIZE(tJac,2)
    DO iEq=nA,nB
       ! compute maximal species contribution for component iC
       xmax = maxval ( abs( tJac(iEq,1:N)) ) 
       DO j=1,N
          ! clean small relative contributions of the line
          x = tJac(iEq,j)
          !WRITE(*,*) "x=", x, xmax
          IF (abs(x)<xeps*xmax) THEN
             !WRITE(*,*) "clean value tJac"
             tJac(iEq,j) = Zero
          END IF
       END DO
    END DO

  END SUBROUTINE Clean_Jac_Section

  !--
  
  SUBROUTINE Clean_Jac_Index(tJac, IdxEq)
    !--------------------------------------------------------
    ! Clean Indexed Equations
    !--------------------------------------------------------
    IMPLICIT NONE
    REAL(dp) :: tJac(:,:)
    INTEGER  :: IdxEq(:)
    !---
    REAL(kind=8) :: x, xmax
    INTEGER  :: iEq,j
    INTEGER  :: N, nEq
    !CALL Info_( "Clean_Jac_Index" )
    N   = SIZE(tJac,2)
    nEq = SIZE(IdxEq)
    
    DO iEq = IdxEq(1),IdxEq(nEq)
       ! compute maximal species contribution for component iC
       xmax = maxval ( abs( tJac(iEq,1:N)) ) 
       DO j=1,N
          ! clean small relative contributions of the line
          x = tJac(iEq,j)
          IF (abs(x)<xeps*xmax) tJac(iEq,j) = Zero
       END DO
    END DO

  END SUBROUTINE Clean_Jac_Index
  
END MODULE M_Clean_Jac
