module M_Matrix_Utils

  use M_Kinds
  implicit none
  private
  real(dp), parameter :: default_epsilon = 1.D-10
  
  public :: Matrix_Clear_Zeros
  public :: Matrix_Norm_L1
  public :: Matrix_Norm_L2
  public :: Matrix_Norm_Linf
  public :: Matrix_Error

contains

!---

subroutine Matrix_Clear_Zeros(tMat, user_epsilon) !NEW_17/10/2007 17:37
  !==========================================================
  ! Purpose : Clear Small Matrix Coefficients (Line By Line)
  !==========================================================
  implicit none
  real(dp) :: tMat(:,:)
  integer :: i, j
  integer :: nl, nk
  real(dp) :: maxline
  real(dp), optional :: user_epsilon
  !---
  real(dp) :: cut_epsilon
  !---
  nl = size(tMat,1)
  nk = size(tMat,2)

  cut_epsilon = default_epsilon
  if ( present(user_epsilon)) cut_epsilon = user_epsilon

  do i= 1,nl
     maxline = maxval ( abs(tMat(i,1:nk)) )
     do j= 1,nk
        if ( abs(tMat(i,j)) < cut_epsilon*maxline ) tMat(i,j) = Zero
     end do
  end do  

end subroutine Matrix_Clear_Zeros


  !---
  
  function Matrix_Norm_L1(tMat) result (norm)
    real(dp) :: tMat(:,:)
    real(dp) :: norm, h
    integer :: Nl, Nk, i, j
    !
    Nl = size(tMat,1)
    Nk = size(tMat,2)
    
    h  = 1.D0/max(Nl,Nk)    ! norm(Id)=1
    
    norm = Zero
    do i = 1,Nl
       do j = 1,Nk
       norm = norm + h*abs(tMat(i,j)) 
       end do
    end do
    
  end function Matrix_Norm_L1

  !---
  
  function Matrix_Norm_L2(tMat) result ( norm)
    real(dp) :: tMat(:,:)
    real(dp) :: norm, h
    integer :: Nl, Nk, i, j
    !
    Nl = size(tMat,1)
    Nk = size(tMat,2)
  
    h = 1.D0/max(Nl,Nk) ! norm(Id)=1

    norm = Zero
    do i = 1,Nl
       do j = 1,Nk
       norm = norm + h*(tMat(i,j)*tMat(i,j))
       end do
    end do
    norm = sqrt(norm)
    
  end function Matrix_Norm_L2

  !---
  
  function Matrix_Norm_Linf(tMat) result (norm)
    real(dp) :: tMat(:,:)
    real(dp) :: norm
    integer :: Nl, Nk, i, j
    !
    Nl = size(tMat,1)
    Nk = size(tMat,2)
    norm = Zero
    do i = 1,Nl
       do j = 1,Nk
          norm = max(norm,abs(tMat(i,j)))
       end do
    end do
    
  end function Matrix_Norm_Linf

  !---

  subroutine Matrix_Error(tMat, Filter, tResult) 
    real(dp), intent(in)  :: tMat(:,:)
    real(dp), intent(in)  :: Filter
    real(dp), intent(out) :: tResult(:,:)
    real(dp) :: value
    integer :: Nl, Nk, i, j
    !
    Nl = size(tMat,1)
    Nk = size(tMat,2)
   
    do i = 1,Nl
       do j = 1,Nk
          value = abs(tMat(i,j))
          if (value > filter) then 
             tResult(i,j) = value
          else
             tResult(i,j) = Zero
          end if
       end do
    end do
    
  end subroutine Matrix_Error


end module M_Matrix_Utils
