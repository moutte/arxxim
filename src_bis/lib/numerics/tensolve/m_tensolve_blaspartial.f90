module M_TenSolve_BlasPartial
  ! Contains only DNRM2, DSWAP, IDAMAX & DGEMV
  ! This very much simplified BLAS module for use by TOMS 768 is by Alan Miller
  ! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj
  ! Latest revision - 19 January 1999
  use M_Kinds
  implicit none
  !
  private
  !
  public:: dnrm2
  public:: idamax
  public:: dgemv
  public:: dswap
  !
contains

function dnrm2 ( n, x, incx) result(fn_val)
  !  Euclidean norm of the n-vector stored in x() with storage increment incx .
  !  if n <= 0 return with result= 0.
  !  if n >= 1 then incx must be >= 1
  !  c.l.lawson, 1978 jan 08
  !  modified to correct failure to update ix, 1/25/92.
  !  modified 3/93 to return if incx <= 0.
  !  This version by Alan.Miller @ vic.cmis.csiro.au
  !  Latest revision - 22 January 1999
  !  four phase method using two built-in constants that are
  !  hopefully applicable to all machines.
  !      cutlo= maximum of  SQRT(u/eps)  over all known machines.
  !      cuthi= minimum of  SQRT(v)      over all known machines.
  !  where
  !      eps= smallest no. such that eps + 1. > 1.
  !      u  = smallest positive no.   (underflow limit)
  !      v  = largest  no.            (overflow  limit)
  !  brief outline of algorithm..
  !  phase 1    scans zero components.
  !  move to phase 2 when a component is nonzero and <= cutlo
  !  move to phase 3 when a component is > cutlo
  !  move to phase 4 when a component is >= cuthi/m
  !  where m= n for x() real and m= 2*n for complex.
  integer, intent(in):: n, incx
  real(dp),intent(in):: x(:)
  real(dp)           :: fn_val
  ! Local variables
  integer :: i, ix, j, next
  real(dp):: cuthi, cutlo, hitest, sum, xmax
  if(n <= 0 .or. incx <= 0) then
    fn_val= zero
    return
  end if
  !_______________________________Set machine-dependent constants
  cutlo= SQRT( TINY(one) / EPSILON(one) )
  cuthi= SQRT( HUGE(one) )
  next= 1
  sum= zero
  i= 1
  ix= 1
  !_______________________________begin main loop
  20 select case (next)
    case (1)
       if( ABS(x(i)) > cutlo) GO TO 85
       next= 2
       xmax= zero
       GO TO 20
    case (2)
       !__________________________phase 1.  sum is zero
       if( x(i) == zero) GO TO 200
       if( ABS(x(i)) > cutlo) GO TO 85
       !__________________________prepare for phase 2.   x(i) is very small.
       next= 3
       GO TO 105
    case (3)
       !__________________________phase 2.  sum is small.
       !__________________________scale to avoid destructive underflow.
       if( ABS(x(i)) > cutlo ) then
         !prepare for phase 3.
         sum= (sum * xmax) * xmax
         GO TO 85
       end if
    case (4)
       GO TO 110
  end select
  !
  !_______________________________common code for phases 2 and 4.
  !_______________________________in phase 4 sum is large.  scale to avoid overflow.
  110 if( ABS(x(i)) <= xmax ) GO TO 115
  sum= one + sum * (xmax / x(i))**2
  xmax= ABS(x(i))
  GO TO 200
  !
  !_______________________________phase 3.  sum is mid-range.  no scaling.
  !_______________________________for real or d.p. set hitest= cuthi/n
  !_______________________________for complex      set hitest= cuthi/(2*n)
  85 hitest= cuthi / real( n, dp )
  do j= ix, n
    if(ABS(x(i)) >= hitest) GO TO 100
    sum= sum + x(i)**2
    i= i + incx
  end do
  fn_val= SQRT( sum )
  return
  !
  !_______________________________prepare for phase 4.
  !_______________________________ABS(x(i)) is very large
  100 ix= j
  next= 4
  sum= (sum / x(i)) / x(i)
  !
  !Set xmax; large if next= 4, small if next= 3
  105 xmax= ABS(x(i))
  115 sum= sum + (x(i)/xmax)**2
  200 ix= ix + 1
  i= i + incx
  if( ix <= n ) GO TO 20
  !_______________________________end of main loop.
  !_______________________________compute square root and adjust for scaling.
  fn_val= xmax *SQRT(sum)
  return
end function dnrm2

subroutine dswap (n, x, incx, y, incy)
  !interchanges two vectors.
  integer, intent(in)       :: n, incx, incy
  real(dp),intent(IN OUT) :: x(:), y(:)
  ! Local variables
  real(dp) :: temp(n)
  if(n <= 0) return
  if(incx == 1 .and. incy == 1) then
    temp=  x(:n)
    x(:n)= y(:n)
    y(:n)= temp
    return
  end if
  temp= x(:n*incx:incx)
  x(:n*incx:incx)= y(:n*incy:incy)
  y(:n*incy:incy)= temp
  return
end subroutine dswap

function idamax(n, x, incx) result(fn_val)
  !     finds the index of element having max. absolute value.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  integer, intent(in):: n, incx
  real(dp),intent(in):: x(:)
  !
  integer:: fn_val
  integer:: imax(1)
  fn_val= 0
  if( n < 1 .or. incx <= 0 ) return
  fn_val= 1
  if(n == 1) return
  if(incx == 1) then; imax= MAXLOC( ABS(x(:n)) )
  else             ; imax= MAXLOC( ABS(x(:n*incx:incx)) )
  end if
  fn_val= imax(1)
  return
end function idamax

subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
  ! ELF90 translation by Alan Miller   31-Aug-1997
  !     .. Scalar Arguments ..
  real(dp),intent(in)           :: alpha, beta
  integer, intent(in)           :: incx, incy, lda, m, n
  character (len=1), intent(in) :: trans
  !     .. Array Arguments ..
  real(dp),intent(in)           :: a(:,:), x(:)
  real(dp),intent(IN OUT)       :: y(:)
  !     ..
  !  Purpose
  !  == == ===
  !  DGEMV  performs one of the matrix-vector operations
  !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  !  where alpha and beta are scalars, x and y are vectors and A is an
  !  m by n matrix.
  !  Parameters
  !  == == == == ==
  !  TRANS  - character*1.
  !           On entry, TRANS specifies the operation to be performed as
  !           follows:
  !              TRANS= 'N' or 'n'   y := alpha*A*x + beta*y.
  !              TRANS= 'T' or 't'   y := alpha*A'*x + beta*y.
  !              TRANS= 'C' or 'c'   y := alpha*A'*x + beta*y.
  !           Unchanged on exit.
  !  M      - integer.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !  N      - integer.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !  ALPHA  - doUBLE PRECISION.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !  A      - doUBLE PRECISION array of dimension ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !  LDA    - integer.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least max( 1, m ).
  !           Unchanged on exit.
  !  X      - doUBLE PRECISION array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS= 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !  INCX   - integer.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !  BETA   - doUBLE PRECISION.
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !  Y      - doUBLE PRECISION array of dimension at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS= 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
  !  INCY   - integer.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !  Level 2 Blas routine.
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !     .. Local Scalars ..
  real(dp):: temp
  integer ::  i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
  !     .. External Functions ..
  ! logical ::  lsame
  ! external           lsame
  !     .. External Subroutines ..
  ! external           erinfo
  !     .. Intrinsic Functions ..
  ! intrinsic          MAX
  !     ..
  !     .. Executable Statements ..
  !     Test the input parameters.
  info= 0
  if (    trans/='N' .and. trans/='n' &
  & .and. trans/='T' .and. trans/='t' &
  & .and. trans/='C' .and. trans/='c') then; info= 1
  elseif(m < 0)             then; info= 2
  elseif(n < 0)             then; info= 3
  elseif(lda < MAX( 1, m )) then; info= 6
  elseif(incx  == 0)         then; info= 8
  elseif(incy  == 0)         then; info= 11
  end if
  if( info /= 0 ) then
    write(*, '(a, i4, a)') ' Error number: ', info, ' in BLAS2 routine DGEMV'
    return
  end if
  !_______________________________________Quick return if possible.
  if(     (m ==0) .or. (n ==0) &
  &  .or. ((alpha ==zero) .and. (beta ==one))) return
  !Set  lenX  and  lenY, the lengths of the vectors x and y, and set
  !up the start points in  X  and  Y.
  if (trans=='N' .or. trans=='n') then
    lenx= n; leny= m
  else
    lenx= m; leny= n
  end if
  if( incx > 0 ) then; kx= 1
  else               ; kx= 1 - ( lenx - 1 )*incx
  end if
  if( incy > 0 ) then; ky= 1
  else               ; ky= 1 - ( leny - 1 )*incy
  end if
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !     First form  y := beta*y.
  if( beta /= one ) then
    if( incy == 1 ) then
      if( beta == zero ) then
        y( :leny )= zero
      else
        y( :leny )= beta*y( :leny )
      end if
    else
      iy= ky
      if( beta == zero ) then
        do i= 1, leny
          y( iy )= zero
          iy     = iy   + incy
        end do
      else
        do i= 1, leny
          y( iy )= beta*y( iy )
          iy     = iy           + incy
        end do
      end if
    end if
  end if
  if( alpha == zero ) return
  if ( trans == 'N' .or. trans == 'n' ) then
  !        Form  y := alpha*A*x + y.
    jx= kx
    if( incy == 1 ) then
      do j= 1, n
        if( x( jx ) /= zero ) then
          temp= alpha*x( jx )
          y( 1:m )= y( 1:m ) + temp*a( 1:m, j )
        end if
        jx= jx + incx
      end do
    else
      do j= 1, n
        if( x( jx ) /= zero ) then
          temp= alpha*x( jx )
          iy  = ky
          do i= 1, m
            y( iy )= y( iy ) + temp*a( i, j )
            iy     = iy      + incy
          end do
        end if
        jx= jx + incx
      end do
    end if
  else
  !        Form  y := alpha*A'*x + y.
    jy= ky
    if( incx == 1 ) then
      do j= 1, n
        temp= dot_product( a(1:m,j), x(1:m) )
        y( jy )= y( jy ) + alpha*temp
        jy     = jy      + incy
      end do
    else
      do j= 1, n
        temp= zero
        ix  = kx
        do i= 1, m
          temp= temp + a( i, j )*x( ix )
          ix  = ix   + incx
        end do
        y( jy )= y( jy ) + alpha*temp
        jy     = jy      + incy
      end do
    end if
  end if
  return
  !     End of DGEMV.
end subroutine dgemv

end module M_TenSolve_BlasPartial
