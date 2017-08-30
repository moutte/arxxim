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

subroutine Simplex_Calc(iError) ! &
!! & tSimplex,IZROV,IPOSV, &
!! & iError,n1,n2)
!--
!-- simplex routine, 
!-- modified from Press et al, Numerical Recipes in Fortran
!--
  use M_Simplex_Vars, only: tSimplex,IZROV,IPOSV
  use M_Numeric_Tools,only: iFirstLoc,iMaxLoc_R,OuterProd_R, Swap_I
  
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
  INIT=    .true.
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

