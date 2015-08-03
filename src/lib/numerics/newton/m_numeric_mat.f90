MODULE M_Numeric_Mat
!.routines from Numerical Recipes"
  USE M_Kinds
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC:: LU_Decomp
  PUBLIC:: LU_BakSub
  PUBLIC:: CheckSingular
  !
CONTAINS

SUBROUTINE CheckSingular( & !"NR"
& A, &       !IN:  a square matrix
& N, &       !IN:  INTEGER, DIMENSION of upper-left block of A 
& Singular)  !OUT: LOGICAL
  USE M_Numeric_Tools,   ONLY: iMaxLoc_R,OuterProd_R,Swap_RV
  REAL(dp),INTENT(IN) :: A(:,:)
  INTEGER, INTENT(IN) :: N
  LOGICAL, INTENT(OUT):: Singular
  !
  REAL(dp),PARAMETER:: TINY=1.0E-18_dp
  REAL(dp):: W(N)  !stores the IMPLICIT scaling of each row.
  REAL(dp):: AA(N,N)
  INTEGER :: J,iMax !,N
  !
  Singular=.FALSE.
  AA(:,:)= A(1:N,1:N)
  !N= SIZE(A,1)
  !
  W=MAXVAL(ABS(AA),DIM=2) !Loop over rows to get the IMPLICIT scaling information
  IF (ANY(W==Zero)) THEN !there is at least one row with ONLY Zeros ...
    Singular= .TRUE.
    RETURN
  ENDIF
  W= One/W !Save the scaling.
  DO J=1,N
    iMax= (J-1) + iMaxLoc_R(W(J:N)*ABS(AA(J:N,J))) !Find the pivot row.
    IF(J /= iMAx) THEN !need to interchange rows ?
      CALL Swap_RV(AA(iMax,:),AA(J,:)) !interchange,
      W(iMax)=W(J) !also interchange the scale factor
    ENDIF
    !! IF (AA(J,J)==Zero) THEN
    IF (ABS(AA(J,J))<Tiny) THEN
      Singular=.TRUE.
      RETURN
    ENDIF
    AA(J+1:N,J)=    AA(J+1:N,J) /AA(J,J) !Reduce remaining submatrix.
    AA(J+1:N,J+1:N)=AA(J+1:N,J+1:N) -OuterProd_R(AA(J+1:N,J),AA(J,J+1:N))
  ENDDO
  RETURN
ENDSUBROUTINE CheckSingular

SUBROUTINE LU_Decomp( & !"NR"
& A,    &    !IN:  a matrix, OUT: LU decompsition of a rowwise permutation of itself
& Indx, &    !OUT: records row permutations effected by the partial pivoting
& D,    &    !OUT: +1/-1 = even/odd number of row interchanges
& Singular)  !OUT
!USEd in combination with LU_BakSub to solve linear equations or invert a matrix
  !
  USE M_Numeric_Tools,   ONLY: iMaxLoc_R,OuterProd_R,Swap_RV
  USE M_Trace,ONLY: Stop_
  REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
  INTEGER, DIMENSION(:),  INTENT(OUT)  :: Indx
  REAL(dp),               INTENT(OUT)  :: D
  LOGICAL,                INTENT(OUT)  :: Singular
  !
  REAL(dp),DIMENSION(SIZE(A,1)) :: W  !stores the IMPLICIT scaling of each row.
  REAL(dp),PARAMETER :: TINY=1.0E-20_dp
  INTEGER:: J,N,IMAX
  !
  Singular=.FALSE.
  N=SIZE(Indx)
  IF (SIZE(A,1)/=N .OR. SIZE(A,2)/=N) CALL Stop_("error on DIMENSIONs in LU_Decomp")
  !
  D=One !No row interchanges yet.
  W=MAXVAL(ABS(A),DIM=2) !Loop over rows to get the IMPLICIT scaling information
  IF (ANY(W==Zero)) THEN !there is at least one row with ONLY Zeros ...
    Singular=.TRUE.
    RETURN
  ENDIF
  W=One/W !Save the scaling.
  DO J=1,N
    IMAX=(J-1)+iMaxLoc_R(W(J:N)*ABS(A(J:N,J))) !Find the pivot row.
    IF(J /= IMAX) THEN !______________Do we need to interchange rows ?
      CALL Swap_RV(A(IMAX,:),A(J,:)) !interchange,
      D=-D !__________________________and change the parity of d.
      W(IMAX)=W(J) !__________________Also interchange the scale factor
    ENDIF
    Indx(J)=IMAX
    !!!§§§ IF (A(J,J)==Zero) A(J,J)=TINY
    IF (A(J,J)==Zero) THEN
      Singular=.TRUE.
      RETURN
    ENDIF
    !If the pivot element is zero the matrix is singular
    !(at least to the precision of the algorithm,
    ! for some applications on singular matrices, 
    ! it is desirable to substitute TINY for zero)
    A(J+1:N,J)=    A(J+1:N,J) /A(J,J) !Reduce remaining submatrix.
    A(J+1:N,J+1:N)=A(J+1:N,J+1:N) -OuterProd_R(A(J+1:N,J),A(J,J+1:N))
  ENDDO
ENDSUBROUTINE LU_Decomp

SUBROUTINE LU_BakSub( & !"NR"
!.solves the system of linear equations A·X=  B
& A,    & !IN,    matrix, LU decomposition, determined by the routine LU_Decomp  
& Indx, & !IN,    permutation vector RETURNed by LU_Decomp
& B)      !INOUT, IN=the right-hand-side vector, OUT=the solution vector
!NOTE FROM "NR"
!A and Indx are not modIFied by this routine
!and can be left in place for successive calls with dIFferent right-hand sides b
!LU_BakSub takes into account the possibility that B will begin with many zero elements, 
!so it is efficient for USE in matrix inversion.
  USE M_Trace,ONLY:Stop_
  REAL(dp),DIMENSION(:,:),INTENT(IN)    :: A
  INTEGER, DIMENSION(:),  INTENT(IN)    :: Indx
  REAL(dp),DIMENSION(:),  INTENT(INOUT) :: B
  !
  INTEGER:: N,I,J,K
  REAL(dp):: Summ
  !
  N=SIZE(Indx)
  IF (SIZE(A,1)/=N .OR. SIZE(A,2)/=N) CALL Stop_("error on DIMENSIONs in LU_BakSub")
  J=0
  DO I=1,N
    K=    Indx(I)
    Summ= B(K)
    B(K)= B(I)
    IF (J/=0) THEN;           Summ=Summ-DOT_PRODUCT(A(I,J:I-1),B(J:I-1))
    ELSEIF (Summ/=Zero) THEN; J=I
    ENDIF
    B(I)=Summ
  ENDDO
  DO I=N,1,-1
    B(I)=(B(I)-DOT_PRODUCT(A(I,I+1:N),B(I+1:N)))/A(I,I)
  ENDDO
ENDSUBROUTINE LU_BakSub

ENDMODULE M_Numeric_Mat

