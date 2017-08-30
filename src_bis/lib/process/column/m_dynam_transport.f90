module m_dynam_transport
  use m_kinds

  implicit none
  
  private
  
  public:: T1D_Adv_Explicit
  public:: T1D_Adv_TVD
  public:: T1D_Adv_Dis_Implicit
  public:: T1D_Dispersion
  
contains

subroutine T1D_Adv_Explicit(vv,dt,dx, Rval, nx, pConcvec, Concvec)
!-----------------------------------------------------------------------
!-- explicit advection with backward difference
!-----------------------------------------------------------------------
  real(dp),intent(in):: vv,dt,dx
  real(dp),intent(in):: Rval
  integer, intent(in):: nx
  real(dp),intent(in):: pConcvec(:)
  real(dp),intent(out):: Concvec(:)
  !
  real(dp):: Cr
  !---------------------------------------------------------------------
  cr = vv * dt / dx
  !
  Concvec(2:nx) = -(cr / Rval) &
  &         * (pConcvec(2:nx) - pConcvec(1:nx-1)) + pConcvec(2:nx)

  return
end subroutine T1D_Adv_Explicit

subroutine T1D_Adv_Dis_Implicit(vv,dt,dx,Disp, Rval, nx, Concvec)
  real(dp),intent(in):: vv,dt,dx,Disp
  real(dp),intent(in):: Rval
  integer, intent(in):: nx
  real(dp),intent(inout):: Concvec(:)
  !
  real(dp),dimension(:),allocatable:: a,b,c,d,v
  real(dp):: alpha,beta
  !
  integer:: i
  !---------------------------------------------------------------------
  alpha = vv * dx /Disp /2.D0
  beta = dx**2 /Disp /dt
  
  allocate(a(nx),b(nx),c(nx),d(nx),v(nx))
  
  ! setting the a, b, c, d values for the tridiagonal matrix
  a(2:nx-1) =  1.D0 + alpha
  b(2:nx-1) = -2.D0 - beta * Rval
  c(2:nx-1) =  1.D0 - alpha

  a(1) = 0.D0
  b(1) = 1.D0
  c(1) = 0.D0

  a(nx) =  alpha *2.D0
  b(nx) = -(beta *Rval + alpha*2.D0)
  c(nx) =  0.D0

  d(1) = Concvec(1)
  d(2:nx) = -beta * Rval * Concvec(2:nx)

  call tridiagonal(nx, a, b, c, d, Concvec)
 
  deallocate(a,b,c,d,v)
 
  return
end subroutine T1D_Adv_Dis_Implicit

subroutine T1D_Dispersion(dt,dx,Disp, Rval, nx, Concvec)
!-----------------------------------------------------------------------
!-- solves the dispersion part of the equation
!-- using the implicit finite difference
!-----------------------------------------------------------------------
  real(dp),intent(in):: dt,dx,Disp
  real(dp),intent(in):: Rval
  integer, intent(in) :: nx
  real(dp),intent(inout):: Concvec(:)
  !---------------------------------------------------------------------
  real(dp),dimension(nx):: a,b,c,d
  real(dp):: beta,omega
  integer:: i
  !---------------------------------------------------------------------
  omega= 1.D0
  beta= Disp * dt /Rval /dx**2

  a(2:nx-1) = -omega * beta
  b(2:nx-1) =  omega * beta *2.0D0 + 1.D0
  c(2:nx-1) = -omega * beta
  d(2:nx-1) =  Concvec(2:nx-1) &
  &       + (1.D0 - omega) * beta &
  &       * (Concvec(3:nx) - 2.D0 * Concvec(2:nx-1) + Concvec(1:nx-2))
  
  a(1) = 0.D0
  b(1) = 1.D0
  c(1) = 0.D0
  d(1) = Concvec(1)
  
  a(nx) = -omega *beta
  b(nx) = 1.D0 + omega *beta
  c(nx) = 0.D0
  d(nx) = ((1.D0 - omega) *beta *(-Concvec(nx) + Concvec(nx-1))) &
  &     + Concvec(nx)
  
  call tridiagonal(nx, a, b, c, d, Concvec)
  
  return
end subroutine T1D_Dispersion

subroutine tridiagonal(n,a,b,c,d,v)
  integer,intent(in):: n
  real(dp),dimension(:),intent(in):: a,b,c,d
  real(dp),dimension(:),intent(out):: v
  
  integer:: i, j
  real(dp),dimension(n):: b1,d1,ff
  
  b1(:) = b(:)
  d1(:) = d(:)
  
  do i= 2,n
    ff(i) = a (i) / b1(i - 1)
    b1(i) = b1(i) - c (i - 1)  *ff(i)
    d1(i) = d1(i) - d1(i - 1) *ff(i)
  end do

  v(n) = d1(n) / b1(n)

  do j= n-1,1,-1
    v(j) = (d1(j) - c(j) * v(j + 1)) / b1(j)
  end do

  return
end subroutine tridiagonal

subroutine T1D_Adv_TVD(vv,dt,dx, Rval,nx,c,cnew)
!-----------------------------------------------------------------------
!-- TVD advection using van-leer flux limiter
!-----------------------------------------------------------------------
  real(dp),intent(in):: vv,dt,dx
  real(dp),intent(in):: Rval
  integer, intent(in):: nx
  real(dp),intent(in):: c(:)
  real(dp),intent(out):: cnew(:)
  
  real(dp):: phiL(nx+1)
  real(dp):: phiR(nx+1)
  real(dp):: dc1,dc2,theta
  real(dp):: tol,cr
  integer:: i
  
  tol = 1.D-6
  cr = vv *dt /dx /Rval
    
  phiL(1:nx+1) = 0.D0
  phiR(1:nx+1) = 0.D0
  
  do i = 3,nx
    dc1 = c(i) - c(i - 1)
    dc2 = c(i - 1) - c(i - 2)
    if (abs(dc1) + abs(dc2) < tol) Then
      theta = 1.D0
    elseif (abs(dc1) < tol .and. abs(dc2) > tol) then
      theta = 100.D0
    else
      theta = dc2 / dc1
    end if
    phiL(i) = (theta + abs(theta)) / (1.D0 + abs(theta))
  end do
  
  cnew(2:nx) = c(2:nx) &
  &          - cr *(c(2:nx) - c(1:nx-1)) &
  &          - 0.5D0 *cr *(1.D0 - cr) &
  &           *( phiL(3:nx+1) *(c(3:nx+1) - c(2:nx  )) &
  &            - phiL(2:nx)   *(c(2:nx  ) - c(1:nx-1)) )
  
  !pause
  return
end subroutine T1D_Adv_TVD

end module m_dynam_transport
