module M_Simplex_Calc
!!!  using Simplex_Calc,
!!!  for a system of (nCn Components, nFs Fases)
!!!  where
!!!    vMolCpn(1:nCpn) is the components' abundance
!!!    vFasGibbs(1:nFs) is the array of the fases' free energy
!!!    tStoikSpl is the stoikiometry table of the Fases in terms of the components
!!!
!!!  call Simplex_Vars_Alloc(nCpn,nFs)
!!!  !
!!!  !first row= Gibbs energy of phases 1:nFs
!!!  tSimplex(0,1:nFs)= vFasGibbs(1:nFs)
!!!  !first column= bulk compos'n  
!!!  tSimplex(1:nCpn,0)= vMolCpn(1:nCpn)       
!!!  !main = stoikiometry matrix              
!!!  tSimplex(1:nCpn,1:nFs)=-transpose(tStoikSpl(1:nFs,1:nCpn))
!!!  !                      
!!!  call Simplex_Calc(iError)
!!!  !
!!!  iError is  0 =  Ok
!!!  iError is  1 = "Unbounded Objective Function"
!!!  iError is -1 = "No Solutions Satisfy Constraints Given"
!
  use M_Kinds
  use M_Trace,only:Stop_
  implicit none
  !
  private
  !
  public:: Simplex_Calc
  !
  logical,parameter:: new=.true.
  !
contains

subroutine Simplex_Calc( &
& M,N,                   &
& tSimplex,              &
& IZROV,IPOSV,           &
& iError) !,n1,n2)
!--
!-- simplex routine, 
!-- modified from Press et al, Numerical Recipes in Fortran
!--
  ! use M_Simplex_Vars, only: tSimplex,IZROV,IPOSV
  !
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R,Swap_I
  
  ! M=nCpn
  ! N=nFs
  
  integer, intent(in)   :: M,N
  real(dp),intent(inout):: tSimplex(0:M+1,0:N)
  integer, intent(out)  :: IZROV(N),IPOSV(M)
  integer, intent(out)  :: iError
  !
  real(dp),allocatable:: a(:,:)
  ! integer :: M,N
  
  iError= 0
  !
  ! M= SIZE(tSimplex,1) -2 !or M= SIZE(IPOSV)
  ! N= SIZE(tSimplex,2) -1 !or N= SIZE(IZROV)
  !
  if(new) then
    call simplx(tSimplex,0,0,M,iError,izrov,iposv)
  else  
    allocate(a(M+2,N+1))
    a(1:M+2,1:N+1)= tSimplex(0:M+1,0:N)
    !
    call simplx(a,0,0,M,iError,izrov,iposv)
    !
    tSimplex(0:M+1,0:N)= a(1:M+2,1:N+1)
    deallocate(a)
  end if
  
  return
end subroutine Simplex_Calc

function Arth_I(First,Increment,dimN) !-> array(1:dimN)
!--
!-- build array(dimN) of integer (First,First+Increment,First+2*Increment,...)
!-- ex: Arth_I(1,1,dimN) gives 1,2,3,4,...dimN
!--
  integer, intent(in):: First,Increment,dimN
  integer, dimension(dimN):: Arth_I
  !
  integer :: K,K2,TEMP
  
  if(dimN>0) Arth_I(1)=First
  !if (dimN<=nPar_Arth) then !integer, parameter :: nPar_Arth=16,NPAR2_ARTH=8
  if (dimN<=16) then
    do K=2,dimN
      Arth_I(K)=Arth_I(K-1)+Increment
    end do
    ! Arth_I(2:dimN)=Arth_I(1:dimN-1)+Increment  
  else
  !do K=2,NPAR2_ARTH; Arth_I(K)=Arth_I(K-1)+Increment; end do
    do K=2,8
      Arth_I(K)=Arth_I(K-1)+Increment
    end do
    TEMP=Increment*8 !NPAR2_ARTH
    K=8 !NPAR2_ARTH
    do
      if(K>=dimN) exit
      K2=   K+K
      Arth_I(K+1:MIN(K2,dimN))= TEMP+Arth_I(1:MIN(K,dimN-K))
      Temp= Temp +Temp
      K=  K2
    enddo
  end if
  
  return
end function Arth_I
  
!Simplex, notes from NumRecp:
!
! linear programming= 
! maximize the linear function 
!  Z=  a(0,1).x(1) + ... + a(0,n).x(n)
! subject to 
!  primary constraints: a(0,1:n)>=0
! and simultaneously subject to 
!  M=  m1 + m2 + m3 additional constraints,
!  m1 of them of the form: a(i,1).x1 + ... + a(i,N).x(N) >= b(i)   //i=1..m1
!  m2 of them of the form: a(j,1).x1 + ... + a(j,N).x(N) <= b(j)   //j=m1+1..m1+m2
!  m3 of them of the form: a(k,1).x1 + ... + a(k,N).x(N) <= b(k)   //k=m1+m2+1..m1+m2+m3

! On output, the tableau a is indexed by two returned arrays of integers
! iposv(j) contains, for j= 1...M, 
! the number i whose original variable x(i) is now represented by row j+1 of a
! These are thus the left-hand variables in the solution
! (The first row of a is of course the z-row.)
! A value i > N indicates that the variable is a y(i) rather than an x(i), x(N+j)==y(j)
! Likewise, izrov(j) contains, for j= 1..N,
! the number i whose original variable x(i) is now a right-hand variable, represented by column j+1 of a. 
! These variables are all zero in the solution. 
! The meaning of i > N is the same as above, 
! except that i >N +m1 +m2 denotes an artificial or slack variable 
! which was used only internally and should now be entirely ignored.

subroutine simplx0(a,m1,m2,m3,icase,izrov,iposv)
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R,Swap_I
  implicit none
  
  real(DP), dimension(:,:), intent(inout) :: a
  integer, intent(in) :: m1,m2,m3
  integer, intent(out) :: icase
  integer, dimension(:), intent(out) :: izrov,iposv
  
  real(DP), parameter :: EPS=1.0e-6_dp
  
  integer :: ip,k,kh,kp,nl1,m,n
  integer, dimension(SIZE(a,2)) :: l1
  integer, dimension(m2) :: l3
  real(DP) :: bmax
  logical:: init
  
  !! m=assert_eq(SIZE(a,1)-2,SIZE(iposv),'simplx: m')
  !! n=assert_eq(SIZE(a,2)-1,SIZE(izrov),'simplx: n')
  m=SIZE(a,1)-2
  n=SIZE(a,2)-1
  
  !! if (m /= m1+m2+m3) call nrerror('simplx: bad input constraint counts')
  if (any(a(2:m+1,1) < 0.0)) then
    !! call nrerror('bad input tableau in simplx')
    icase= -2
    return
  end if
  
  nl1=n
  l1(1:n)=Arth_I(1,1,n)
  izrov(:)=l1(1:n)
  iposv(:)=n+Arth_I(1,1,m)
  
  init=.true.
  
  phase1: do
    if (init) then
      if (m2+m3 == 0) exit phase1
      init=.false.
      l3(1:m2)=1
      a(m+2,1:n+1)=-sum(a(m1+2:m+1,1:n+1),dim=1)
    end if
    if (nl1 > 0) then
      kp=l1(imaxloc_R(a(m+2,l1(1:nl1)+1)))
      bmax=a(m+2,kp+1)
    else
      bmax=0.0
    end if
    
    phase1a: do
      if (bmax <= EPS .and. a(m+2,1) < -EPS) then
        icase=-1
        return
      else if (bmax <= EPS .and. a(m+2,1) <= EPS) then
        do ip=m1+m2+1,m
          if (iposv(ip) == ip+n) then
            if (nl1 > 0) then
              kp=l1(imaxloc_R(ABS(a(ip+1,l1(1:nl1)+1))))
              bmax=a(ip+1,kp+1)
            else
              bmax=0.0
            end if
            if (bmax > EPS) exit phase1a
          end if
        end do
        where (spread(l3(1:m2),2,n+1) == 1) &
          a(m1+2:m1+m2+1,1:n+1)=-a(m1+2:m1+m2+1,1:n+1)
        exit phase1
      end if
      call simp1
      if (ip == 0) then
        icase=-1
        return
      end if
      exit phase1a
    end do phase1a
    
    call simp2(m+1,n)
    
    if (iposv(ip) >= n+m1+m2+1) then
      k=ifirstloc(l1(1:nl1) == kp)
      nl1=nl1-1
      l1(k:nl1)=l1(k+1:nl1+1)
    else
      kh=iposv(ip)-m1-n
      if (kh >= 1) then
        if (l3(kh) /= 0) then
          l3(kh)=0
          a(m+2,kp+1)=a(m+2,kp+1)+1.0_dp
          a(1:m+2,kp+1)=-a(1:m+2,kp+1)
        end if
      end if
    end if
    
    call Swap_I(izrov(kp),iposv(ip))
  
  end do phase1
  
  phase2: do
    if (nl1 > 0) then
      kp=l1(imaxloc_R(a(1,l1(1:nl1)+1)))
      bmax=a(1,kp+1)
    else
      bmax=0.0
    end if
    if (bmax <= EPS) then
      icase=0
      return
    end if
    call simp1
    if (ip == 0) then
      icase=1
      return
    end if
    call simp2(m,n)
    call Swap_I(izrov(kp),iposv(ip))
  end do phase2
  
  contains

  subroutine simp1
    integer :: i,k
    real(DP) :: q,q0,q1,qp
    
    ip=0
    i=ifirstloc(a(2:m+1,kp+1) < -EPS)
    if (i > m) return
    
    q1=-a(i+1,1)/a(i+1,kp+1)
    ip=i
    do i=ip+1,m
      if (a(i+1,kp+1) < -EPS) then
        q=-a(i+1,1)/a(i+1,kp+1)
        if (q < q1) then
          ip=i
          q1=q
        else if (q == q1) then
          do k=1,n
            qp=-a(ip+1,k+1)/a(ip+1,kp+1)
            q0=-a(i+1,k+1)/a(i+1,kp+1)
            if (q0 /= qp) exit
          end do
          if (q0 < qp) ip=i
        end if
      end if
    end do
    
  end subroutine simp1
  
  subroutine simp2(i1,k1)
    integer, intent(in) :: i1,k1
    
    integer :: ip1,kp1
    real(DP) :: piv
    integer, dimension(k1) :: icol
    integer, dimension(i1) :: irow
    integer, dimension(MAX(i1,k1)+1) :: itmp
    
    ip1=ip+1
    kp1=kp+1
    
    piv=1.0_dp/a(ip1,kp1)
    
    itmp(1:k1+1)=Arth_I(1,1,k1+1)
    icol=pack(itmp(1:k1+1),itmp(1:k1+1) /= kp1)
    itmp(1:i1+1)=Arth_I(1,1,i1+1)
    irow=pack(itmp(1:i1+1),itmp(1:i1+1) /= ip1)
    
    a(irow,kp1)=a(irow,kp1)*piv
    a(irow,icol)=a(irow,icol)-OuterProd_R(a(irow,kp1),a(ip1,icol))
    a(ip1,icol)=-a(ip1,icol)*piv
    a(ip1,kp1)=piv
    
  end subroutine simp2

end subroutine simplx0

subroutine simplx(a,m1,m2,m3,icase,izrov,iposv)
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R,Swap_I
  implicit none
  !---------------------------------------------------------------------
  real(DP),intent(inout):: a(0:,0:)
  integer, intent(in)   :: m1,m2,m3
  integer, intent(out)  :: icase
  integer, intent(out)  :: izrov(:),iposv(:)
  !---------------------------------------------------------------------  
  real(DP), parameter :: EPS=1.0e-6_dp
  !  
  integer :: ip,k,kh,kp,nl1,m,n
  integer :: l1(size(a,2))
  integer :: l3(m2)
  real(dp):: bmax
  logical :: init
  !---------------------------------------------------------------------  
  !! m=assert_eq(SIZE(a,1)-2,SIZE(iposv),'simplx: m')
  !! n=assert_eq(SIZE(a,2)-1,SIZE(izrov),'simplx: n')
  m=SIZE(a,1)-2
  n=SIZE(a,2)-1
  !! if (m /= m1+m2+m3) call nrerror('simplx: bad input constraint counts')
  if (any(a(1:m,0) < 0.0D0)) then
    !! call nrerror('bad input tableau in simplx')
    icase= -2
    return
  end if
  !
  nl1=n
  l1(1:n)=Arth_I(1,1,n)
  izrov(:)=l1(1:n)
  iposv(:)=n+Arth_I(1,1,m)
  !
  init=.true.
  !
  phase1: do
    if (init) then
      if (m2+m3 == 0) exit phase1
      init=.false.
      l3(1:m2)=1
      a(m+1,0:n)=-sum(a(m1+1:m,0:n),dim=1)
    end if
    if (nl1 > 0) then
      kp=l1(imaxloc_R(a(m+1,l1(1:nl1))))
      bmax=a(m+1,kp)
    else
      bmax=0.0D0
    end if
    
    phase1a: do

      if (bmax <= EPS .and. a(m+1,0) < -EPS) then
        icase=-1
        return
      else if (bmax <= EPS .and. a(m+1,0) <= EPS) then
        do ip=m1+m2+1,m
          if (iposv(ip) == ip+n) then
            if (nl1 > 0) then
              kp=l1(imaxloc_R(ABS(a(ip,l1(1:nl1)))))
              bmax=a(ip,kp)
            else
              bmax=0.0
            end if
            if (bmax > EPS) exit phase1a
          end if
        end do
        where (spread(l3(1:m2),2,n+1) == 1) &
          a(m1+1:m1+m2,0:n)=-a(m1+1:m1+m2,0:n)
        exit phase1
      end if
      call simp1
      if (ip == 0) then
        icase=-1
        return
      end if
      exit phase1a
    end do phase1a
    
    call simp2(m+1,n)
    
    if (iposv(ip) >= n+m1+m2+1) then
      k=ifirstloc(l1(1:nl1) == kp)
      nl1=nl1-1
      l1(k:nl1)=l1(k+1:nl1+1)
    else
      kh=iposv(ip)-m1-n
      if (kh >= 1) then
        if (l3(kh) /= 0) then
          l3(kh)=0
          a(m+1,kp)=a(m+1,kp)+1.0_dp
          a(0:m+1,kp)=-a(0:m+1,kp)
        end if
      end if
    end if
    
    call Swap_I(izrov(kp),iposv(ip))
  
  end do phase1
  
  phase2: do
    !
    if (nl1 > 0) then
      kp=l1(imaxloc_R(a(0,l1(1:nl1))))
      bmax=a(0,kp)
    else
      bmax=0.0
    end if
    if (bmax <= EPS) then
      icase=0
      return
    end if
    call simp1
    if (ip == 0) then
      icase=1
      return
    end if
    call simp2(m,n)
    call Swap_I(izrov(kp),iposv(ip))
    !
  end do phase2
  
  contains

  subroutine simp1
    integer :: i,k
    real(DP):: q,q0,q1,qp
    !
    ip=0
    i=ifirstloc(a(1:m,kp) < -EPS)
    if (i > m) return
    !
    q1=-a(i,0)/a(i,kp)
    ip=i
    do i=ip+1,m
      if (a(i,kp) < -EPS) then
        q=-a(i,0)/a(i,kp)
        if (q < q1) then
          ip=i
          q1=q
        else if (q == q1) then
          do k=1,n
            qp=-a(ip,k)/a(ip,kp)
            q0=-a(i,k) /a(i,kp)
            if (q0 /= qp) exit
          end do
          if (q0 < qp) ip=i
        end if
      end if
    end do
    return
  end subroutine simp1
  
  subroutine simp2(i1,k1)
    integer, intent(in) :: i1,k1
    !
    integer :: ip1,kp1
    real(DP) :: piv
    integer, dimension(k1) :: icol
    integer, dimension(i1) :: irow
    integer, dimension(MAX(i1,k1)+1) :: itmp
    !
    ip1=ip+1
    kp1=kp+1
    !
    piv=1.0_dp/a(ip,kp)
    !
    itmp(1:k1+1)=Arth_I(1,1,k1+1)
    icol=pack(itmp(1:k1+1),itmp(1:k1+1) /= kp1)
    itmp(1:i1+1)=Arth_I(1,1,i1+1)
    irow=pack(itmp(1:i1+1),itmp(1:i1+1) /= ip1)
    !
    irow=irow-1 !JM
    icol=icol-1 !JM
    !
    a(irow,kp)=   a(irow,kp)*piv
    a(irow,icol)= a(irow,icol)-OuterProd_R(a(irow,kp),a(ip,icol))
    a(ip,icol)=  -a(ip,icol)*piv
    a(ip,kp)=     piv
    !
    return
  end subroutine simp2

end subroutine simplx

end module M_Simplex_Calc
