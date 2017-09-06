module M_Simplex_Calc
  !
  use M_Kinds
  use M_Trace,only:Stop_
  !
  implicit none
  !
  private
  !
  public:: Simplex_Calc
  public:: Simplex_Calc_2
  !
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
  !!!  tSimplex(1:nCpn,1:nFs)=-TRANSPOSE(tStoikSpl(1:nFs,1:nCpn))
  !!!  !                      
  !!!  call Simplex_Calc(iError)
  !!!  !
  !!!  iError is  0 =  Ok
  !!!  iError is  1 = "Unbounded Objective Function"
  !!!  iError is -1 = "No Solutions Satisfy Constraints Given"
  !
contains

subroutine Simplex_Calc(iError)
!--
!-- simplex routine, 
!-- modified from Press et al, Numerical Recipes in Fortran
!--
  use M_Simplex_Vars, only: tSimplex,IZROV,IPOSV
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R,Swap_I
  
  ! M=nCpn
  ! N=nFs
  
  !! real(dp),intent(inout):: tSimplex(:,:)
  !! integer, intent(out)  :: IZROV(:),IPOSV(:)
  integer, intent(out)  :: iError !,n1,n2
  !
  integer:: M,N
  integer:: n1,n2
  integer,dimension(:),allocatable:: L1,L2
  !
  real(dp),parameter:: EPS= 1.0E-12_dp  !1.0E-6_dp 
  !
  integer :: NL1,NL2 !IP,K,KH,KP,
  integer :: IP,KP,K,KH !,NL1,NL2,M,N
  real(dp):: bMax
  logical :: INIT,PROCEED
  
  iError= 0
  !
  M= size(tSimplex,1) -2 !or M= size(IPOSV)
  N= size(tSimplex,2) -1 !or N= size(IZROV)
  !
  allocate(L1(1:N+1)) !size(tSimplex,2)
  allocate(L2(1:M+2)) !size(tSimplex,1)
  !
  if (ANY(tSimplex(1:M,0) < Zero)) then
    !! call Stop_('BAD INPUT TABLEAU IN SIMPLX')
    iError= -2
    return
  end if
  !
  NL1=      N
  L1(1:N)=  Arth_I(1,1,N)
  IZROV(:)= L1(1:N)
  !
  NL2=      M
  L2(1:M)=  Arth_I(1,1,M)
  iPosV(:)= N+L2(1:M)
  !
  INIT= .true.
  n1= 0
  n2= 0
  
  DoMain: do
    !------------------------------------------------------------phase 1
    Phase1: do
      
      if (INIT) then
        INIT=.false.
        PROCEED=.false.
        tSimplex(M+1,0:N)= -SUM(tSimplex(1:M,0:N),DIM=1)
      end if
      !
      KP=L1(iMaxLoc_R(tSimplex(M+1,L1(1:NL1))))
      bMax=tSimplex(M+1,KP)
      
      !---------------------------------------------------------phase 1A
      Phase1A: do
        
        if (bMax<=EPS.and.tSimplex(M+1,0)<-EPS) then
          
          iError=-1
          return !------------------------------------------------return
        
        else if (bMax <= EPS .and. tSimplex(M+1,0) <= EPS) then
        
          !M12=1
          do IP=1,M
            if (iPosV(IP)==IP+N) then
              KP=   L1(iMaxLoc_R(ABS(tSimplex(IP,L1(1:NL1)))))
              bMax= tSimplex(IP,KP)
              if (bMax>Zero) exit Phase1A
            end if
          end do
          PROCEED=.true.
          
          exit Phase1
          
        end if
        
        call SIMP1(IP) ; n1= n1 +1
        
        if (IP==0) then
          iError=-1
          return !------------------------------------------------return
        end if
        
        exit Phase1A
        
      enddo Phase1A
      !--------------------------------------------------------/phase 1A
      
      !---------------------------------------------------------phase 1B
      Phase1B: do
        !
        call SIMP2(M+1,N) ; n2= n2 +1
        if (iPosV(IP) >= N+1) then
          K=iFirstLoc(L1(1:NL1)==  KP)
          NL1=NL1-1
          L1(K:NL1)=L1(K+1:NL1+1)
        else
          if (iPosV(IP) < N+1) exit Phase1B
          KH=iPosV(IP)-N
          !if (L3(KH)==0) exit Phase1B
          !L3(KH)=0
        end if
        tSimplex(M+1,KP)=   tSimplex(M+1,KP)  +One
        tSimplex(0:M+1,KP)=-tSimplex(0:M+1,KP)
        exit Phase1B
      enddo Phase1B
      !--------------------------------------------------------/phase 1B
      
      call Swap_I(IZROV(KP),iPosV(IP))
      
      if (PROCEED) exit Phase1
      
    enddo Phase1
    !-----------------------------------------------------------/phase 1
    !
    !------------------------------------------------------------phase 2
    PHASE2: do
      !
      KP=L1(iMaxLoc_R(tSimplex(0,L1(1:NL1))))
      bMax=tSimplex(0,KP)
      
      if (bMax<=Zero) then
        iError=0
        return !--------------------------------------------------return
      end if
      
      call SIMP1(IP)
      n1= n1 +1
      if (IP==0) then
        iError=1
        return !--------------------------------------------------return
      end if
      
      call SIMP2(M,N) ; n2= n2 +1
      call Swap_I(IZROV(KP),iPosV(IP))
      
      if (.not. PROCEED) cycle DoMain
      
    end do PHASE2
    !-----------------------------------------------------------/phase 2
    
    exit DoMain
    
  enddo DoMain
  !
  deallocate(L1,L2) !,L3)
  !
contains
  !
  subroutine SIMP1(IP_)
    use M_Numeric_Tools, only: iFirstLoc
    
    integer, intent(out)::IP_
    
    integer :: I,J,K
    real(dp):: Q,Q0,Q1,QP
    
    IP_=0
    I=iFirstLoc(tSimplex(L2(1:NL2),KP) < -EPS)
    if (I > NL2) return
    
    Q1=-tSimplex(L2(I),0)/tSimplex(L2(I),KP)
    IP_=L2(I)
    do I=I+1,NL2
      
      J=L2(I)
      if (tSimplex(J,KP) < -EPS) then
        Q= -tSimplex(J,0) /tSimplex(J,KP)
        if (Q < Q1) then
          IP_= J
          Q1=  Q
        elseif (Q==Q1) then
          do K=1,N
            QP= -tSimplex(IP_,K) /tSimplex(IP_,KP)
            Q0= -tSimplex(J,  K) /tSimplex(J,  KP)
            if (Q0/=QP) exit
          enddo
          if (Q0<QP) IP_=J
        end if
      end if
      
    enddo
    
    return
  end subroutine SIMP1

  subroutine SIMP2(I1,K1)
    integer, intent(in):: I1,K1
    !
    integer            :: IP1,KP1
    real(dp)           :: PIV
    integer, dimension(K1) :: iCol
    integer, dimension(I1) :: iRow
    integer, dimension(1:MAX(I1,K1)+1):: iTmp
    
    IP1= IP+1
    KP1= KP+1
    !
    PIV= One /tSimplex(IP1-1,KP1-1)
    !
    iTmp(1:K1+1)= Arth_I(1,1,K1+1)
    iCol= PACK(iTmp(1:K1+1),iTmp(1:K1+1) /= KP1)
    !
    iTmp(1:I1+1)=  Arth_I(1,1,I1+1)
    iRow= PACK(iTmp(1:I1+1),iTmp(1:I1+1) /= IP1)
    !
    tSimplex(iRow-1,KP1-1)=  tSimplex(iRow-1,KP1-1)*PIV
    tSimplex(iRow-1,iCol-1)= tSimplex(iRow-1,iCol-1) &
    &        - OuterProd_R(tSimplex(iRow-1,KP1-1),tSimplex(IP1-1,iCol-1))
    tSimplex(IP1-1,iCol-1)= -tSimplex(IP1-1,iCol-1)*PIV
    tSimplex(IP1-1,KP1-1)=   PIV
    !
    return
  end subroutine SIMP2

end subroutine Simplex_Calc

subroutine Simplex_Calc_2( &
& M,N,                   &
& tSimplex,              &
& IZROV,IPOSV,           &
& iError) !,n1,n2)
!--
!-- simplex routine, 
!-- modified from Press et al, Numerical Recipes in Fortran
!--
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R,Swap_I
  ! M=nCpn
  ! N=nFs
  integer, intent(in)   :: M,N
  real(dp),intent(inout):: tSimplex(0:M+1,0:N)
  integer, intent(out)  :: IZROV(N),IPOSV(M)
  integer, intent(out)  :: iError
  !
  iError= 0
  !
  call simplx(tSimplex,0,0,M,iError,izrov,iposv)
  !
  return
end subroutine Simplex_Calc_2

subroutine simplx(a,m1,m2,m3,icase,izrov,iposv)
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R,Swap_I
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(inout):: a(0:,0:)
  integer, intent(in)   :: m1,m2,m3
  integer, intent(out)  :: icase
  integer, intent(out)  :: izrov(:),iposv(:)
  !---------------------------------------------------------------------  
  real(dp), parameter :: EPS=1.0e-6_dp
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
  phase1: do !---------------------------------------------------phase 1
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
    
    phase1a: do !-----------------------------------------------phase 1A 
      if (bmax <= EPS .and. a(m+1,0) < -EPS) then
        icase=-1
        return
      else if (bmax <= EPS .and. a(m+1,0) <= EPS) then
        do ip=m1+m2+1,m
          if (iposv(ip) == ip+n) then
            if (nl1 > 0) then
              kp=l1(imaxloc_R(abs(a(ip,l1(1:nl1)))))
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
    end do phase1a !-------------------------------------------/phase 1A
    
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
  
  end do phase1 !-----------------------------------------------/phase 1
  
  phase2: do !---------------------------------------------------phase 2
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
  end do phase2 !-----------------------------------------------/phase 2
  
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
    integer, dimension(max(i1,k1)+1) :: itmp
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

function Arth_I(First,Increment,dimN) !-> array(1:dimN)
!--
!-- build array(dimN) of integer (First,First+Increment,First+2*Increment,...)
!-- ex: Arth_I(1,1,dimN) gives 1,2,3,4,...dimN
!--
  integer, intent(in)  :: First,Increment,dimN
  integer, dimension(dimN):: Arth_I
  integer :: K,K2,TEMP
  !
  if(dimN>0) Arth_I(1)=First
  !
  !if (dimN<=nPar_Arth) then !integer, parameter :: nPar_Arth=16,NPAR2_ARTH=8
  !
  if (dimN<=16) then
    do K=2,dimN; Arth_I(K)=Arth_I(K-1)+Increment; end do
  else
  !do K=2,NPAR2_ARTH; Arth_I(K)=Arth_I(K-1)+Increment; end do
    do K=2,8; Arth_I(K)=Arth_I(K-1)+Increment; end do
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
  
end module M_Simplex_Calc

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

