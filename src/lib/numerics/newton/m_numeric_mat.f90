module M_Numeric_Mat
!.routines from Numerical Recipes"
  use M_Kinds
  implicit none
  private
  !
  public:: LU_Decomp
  public:: LU_BakSub
  public:: CheckSingular
  !
contains

subroutine CheckSingular( & !"NR"
& A, &       !IN:  a square matrix
& N, &       !IN:  integer, dimension of upper-left block of A 
& Singular)  !OUT: logical
  use M_Numeric_Tools,   only: iMaxLoc_R,OuterProd_R,Swap_RV
  real(dp),intent(in) :: A(:,:)
  integer, intent(in) :: N
  logical, intent(out):: Singular
  !
  real(dp),parameter:: TINY=1.0E-18_dp
  real(dp):: W(N)  !stores the implicit scaling of each row.
  real(dp):: AA(N,N)
  integer :: J,iMax !,N
  !
  Singular=.false.
  AA(:,:)= A(1:N,1:N)
  !N= size(A,1)
  !
  W=maxval(abs(AA),DIM=2) !Loop over rows to get the implicit scaling information
  if (any(W==Zero)) then !there is at least one row with only Zeros ...
    Singular= .true.
    return
  end if
  W= One/W !Save the scaling.
  do J=1,N
    iMax= (J-1) + iMaxLoc_R(W(J:N)*abs(AA(J:N,J))) !Find the pivot row.
    if(J /= iMAx) then !need to interchange rows ?
      call Swap_RV(AA(iMax,:),AA(J,:)) !interchange,
      W(iMax)=W(J) !also interchange the scale factor
    end if
    !! if (AA(J,J)==Zero) then
    if (abs(AA(J,J))<Tiny) then
      Singular=.true.
      return
    end if
    AA(J+1:N,J)=    AA(J+1:N,J) /AA(J,J) !Reduce remaining submatrix.
    AA(J+1:N,J+1:N)=AA(J+1:N,J+1:N) -OuterProd_R(AA(J+1:N,J),AA(J,J+1:N))
  end do
  return
end subroutine CheckSingular

subroutine LU_Decomp( & !"NR"
& A,    &    !IN:  a matrix, OUT: LU decompsition of a rowwise permutation of itself
& Indx, &    !OUT: records row permutations effected by the partial pivoting
& D,    &    !OUT: +1/-1 = even/odd number of row interchanges
& Singular)  !OUT
!used in combination with LU_BakSub to solve linear equations or invert a matrix
  !
  use M_Numeric_Tools,   only: iMaxLoc_R,OuterProd_R,Swap_RV
  use M_Trace,only: Stop_
  real(dp),dimension(:,:),intent(inout):: A
  integer, dimension(:),  intent(out)  :: Indx
  real(dp),               intent(out)  :: D
  logical,                intent(out)  :: Singular
  !
  real(dp),dimension(size(A,1)) :: W  !stores the implicit scaling of each row.
  real(dp),parameter :: TINY=1.0E-20_dp
  integer:: J,N,IMAX
  !
  Singular=.false.
  N=size(Indx)
  if (size(A,1)/=N .or. size(A,2)/=N) call Stop_("error on dimensions in LU_Decomp")
  !
  D=One !No row interchanges yet.
  W=maxval(abs(A),DIM=2) !Loop over rows to get the implicit scaling information
  if (any(W==Zero)) then !there is at least one row with only Zeros ...
    Singular=.true.
    return
  end if
  W=One/W !Save the scaling.
  do J=1,N
    IMAX=(J-1)+iMaxLoc_R(W(J:N)*abs(A(J:N,J))) !Find the pivot row.
    if(J /= IMAX) then !______________Do we need to interchange rows ?
      call Swap_RV(A(IMAX,:),A(J,:)) !interchange,
      D=-D !__________________________and change the parity of d.
      W(IMAX)=W(J) !__________________Also interchange the scale factor
    end if
    Indx(J)=IMAX
    !!!§§§ if (A(J,J)==Zero) A(J,J)=TINY
    if (A(J,J)==Zero) then
      Singular=.true.
      return
    end if
    !If the pivot element is zero the matrix is singular
    !(at least to the precision of the algorithm,
    ! for some applications on singular matrices, 
    ! it is desirable to substitute TINY for zero)
    A(J+1:N,J)=    A(J+1:N,J) /A(J,J) !Reduce remaining submatrix.
    A(J+1:N,J+1:N)=A(J+1:N,J+1:N) -OuterProd_R(A(J+1:N,J),A(J,J+1:N))
  end do
end subroutine LU_Decomp

subroutine LU_BakSub( & !"NR"
!.solves the system of linear equations A·X=  B
& A,    & !IN,    matrix, LU decomposition, determined by the routine LU_Decomp  
& Indx, & !IN,    permutation vector returned by LU_Decomp
& B)      !INOUT, IN=the right-hand-side vector, OUT=the solution vector
!NOTE FROM "NR"
!A and Indx are not MODIFIED by this routine
!and can be left in place for successive calls with different right-hand sides b
!LU_BakSub takes into account the possibility that B will begin with many zero elements, 
!so it is efficient for use in matrix inversion.
  use M_Trace,only:Stop_
  real(dp),dimension(:,:),intent(in)    :: A
  integer, dimension(:),  intent(in)    :: Indx
  real(dp),dimension(:),  intent(inout) :: B
  !
  integer:: N,I,J,K
  real(dp):: Summ
  !
  N=size(Indx)
  if (size(A,1)/=N .or. size(A,2)/=N) call Stop_("error on dimensions in LU_BakSub")
  J=0
  do I=1,N
    K=    Indx(I)
    Summ= B(K)
    B(K)= B(I)
    if (J/=0) then;           Summ=Summ-dot_product(A(I,J:I-1),B(J:I-1))
    elseif (Summ/=Zero) then; J=I
    end if
    B(I)=Summ
  end do
  do I=N,1,-1
    B(I)=(B(I)-dot_product(A(I,I+1:N),B(I+1:N)))/A(I,I)
  end do
end subroutine LU_BakSub

end module M_Numeric_Mat

