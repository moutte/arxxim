module M_Numeric_Tools 

  !---------------------------------------------------------
  ! "NR" Tools from Numerical Recipes"
  !---------------------------------------------------------

  use M_Kinds
  use M_Trace, only:Stop_
  implicit none
 
  private
  
  !---
  public:: Interpol
  public:: linear_interp_2d
  !---
  public:: vAbs
  !---
  public:: Swap_RV
  public:: Swap_I
  !---
  public:: OuterProd_R
  public:: OuterAnd
  !---
  public:: iMaxLoc_R
  public:: iMinLoc_R
  public:: iFirstLoc
  !---
  public:: Unit_Matrix
  public:: Sort
  public:: CalcPermut
  !---
  public:: Jacobian_Numeric
  public:: Jacobian_Numeric_Bis
  !---

  integer,public:: fNewt_I= 0 !counter used for DebNewt
  integer,public:: fNewtF=  0 !"trace" file for Newton, values of the unknown vector at each step, opened if DebNewt
  integer,public:: fNewtR=  0 !"trace" file for Newton, values of the residue at each step, opened if DebNewt 
  integer,public:: fNewtJ=  0 !
  !
  !integer,parameter::NPAR_ARTH=16,NPAR2_ARTH=8 !for ARTH_I

contains


!---------------------------------------------LINEAR TABLE INTERPOLATION


  real(dp) function Interpol(T,iota,N,vT,vP)
  !---------------------------------------------------------------------
  ! provisional, basic linear interpolation of P from a table (vT,vP)
  ! assumes the vT-series is strictly increasing from 1 to N
  !---------------------------------------------------------------------
    implicit none
    real(dp), intent(in):: T,iota
    integer,  intent(in):: N
    real(dp), intent(in):: vT(:),vP(:)
    !
    real(dp):: P,K1,K2
    integer :: I
    !
    if(T<vT(1) .or. T>vT(N)) call Stop_("T outside T-series")
    !
    I=1; do while(T>vT(I)); I=I+1; end do !-> find upper bound
    !
    if(abs(T-vT(I))<=Iota) then
      P= vP(I)
    else
      K1= (T-vT(I)  ) /(vT(I-1) - vT(I)  )
      K2= (T-vT(I-1)) /(vT(I)   - vT(I-1))
      P=  K1 *vP(I-1) + K2 *vP(I)
    end if
    Interpol= P
    return
  end function Interpol


!-----------------------------------------------------------VECTOR UTILS


  real(dp) function vAbs(V) 
    !----------------------
    ! Norm of Vector V
    !----------------------
    implicit none
    real(dp),dimension(:),intent(in) :: V
    vAbs=SQRT(dot_product(V,V))
  end function vAbs

  !---

  subroutine Swap_RV(A,B) 
    !-----------------------
    ! Swap RealVector
    !-----------------------
    implicit none
    real(dp),dimension(:),intent(inout):: A,B
    real(dp),dimension(size(A))        :: DUM
    DUM=A; A=B; B=DUM
  end subroutine Swap_RV

  !---

  function iMaxLoc_R(ARR) 
    !-------------------------------
    ! Localize maximum value
    !-------------------------------
    implicit none
    real(dp),dimension(:), intent(in) :: ARR
    integer               :: iMaxLoc_R
    integer, dimension(1) :: IMAX
    IMAX=MAXLOC(ARR(:))
    iMaxLoc_R=IMAX(1)
  end function iMaxLoc_R

  !---

  function iMinLoc_R(ARR) 
    !-------------------------------
    ! Localize minimum value
    !-------------------------------
    implicit none
    real(dp),dimension(:), intent(in) :: ARR
    integer               :: iMinLoc_R
    integer, dimension(1) :: IMin
    IMin=MINLOC(ARR(:))
    iMinLoc_R=IMin(1)
  end function iMinLoc_R

  !---

  function iFirstLoc(MASK) 
    !-------------------------------
    ! Localize first occurence
    !-------------------------------
    implicit none
    integer                        :: iFirstLoc
    logical,dimension(:),intent(in):: MASK
    !
    integer,dimension(1):: LOC
    !
    LOC=MAXLOC(MERGE(1,0,MASK))
    iFirstLoc=LOC(1)
    if (.not. MASK(ifIRSTLOC)) ifIRSTLOC=size(MASK)+1
  end function ifIRSTLOC

  !---

  function OuterProd_R(A,B) 
    !-------------------------------
    ! Outer product of vectors
    !-------------------------------
    implicit none
    real(dp),dimension(:), intent(in) :: A,B
    real(dp),dimension(size(A),size(B)) :: OuterProd_R
    !
    OuterProd_R= SPread(A,DIM=2,NCOPIES=size(B)) &
         &          * SPread(B,DIM=1,NCOPIES=size(A))
  end function OuterProd_R

  !---

  function OuterAND(A,B) 
    !-------------------------------
    ! Outer operator and 
    !-------------------------------
    implicit none
    logical, dimension(:), intent(in)   :: A,B
    logical, dimension(size(A),size(B)) :: OuterAND
    OuterAND= SPread(A,DIM=2,NCOPIES=size(B)) &
         &   .and. SPread(B,DIM=1,NCOPIES=size(A))
  end function OuterAND

  !---

  subroutine Unit_Matrix(MAT) 
    !------------------------------------------------------------
    ! Square real matrix Identity, dimension=MIN(nRows,nColumns)
    !------------------------------------------------------------
    implicit none
    real(dp),dimension(:,:),intent(out) :: MAT
    integer:: I,N
    N=MIN(size(MAT,1),size(MAT,2))
    MAT(:,:)=Zero
    do I=1,N
       MAT(I,I)=One
    end do
  end subroutine UNIT_MATRIX


!----------------------------------------------------------------SORTING

  subroutine SWAP_R(A,B) 
    !-----------------------
    ! Swap Real
    !-----------------------
    implicit none
    real(dp),intent(inout) :: A,B
    real(dp) :: DUM
    DUM=A; A=B; B=DUM
  end subroutine SWAP_R

  !---

  subroutine SWAP_I(A,B)
    !-----------------------
    ! Swap Integer
    !-----------------------
    implicit none
    integer,intent(inout) :: A,B
    integer:: DUM
    DUM=A; A=B; B=DUM
  end subroutine SWAP_I

  !---

  subroutine SWAP_RMASKED(A,B,MASK)
    !-----------------------
    ! Swap Real with Mask
    !-----------------------
    implicit none
    real(dp),intent(inout) :: A,B
    logical,  intent(in)   :: MASK
    real(dp):: SWP
    if (MASK) then
       SWP=A; A=B; B=SWP
    end if
  end subroutine SWAP_RMASKED

  !---

  subroutine Sort(arr) 
    !-------------------------------------------------------------
    ! NR, in ModSort Quicksort
    !------
    ! sorts an array arr into ascending numerical order using 
    ! the Quicksort algorithm. 
    !------
    ! arr is replaced on output by its sorted rearrangement.
    ! Parameters:
    ! NN=     size of subarrays sorted by straight insertion
    ! NSTACK= required auxiliary storage.
    !-------------------------------------------------------------
    implicit none
    real(dp),dimension(:), intent(inout) :: arr
    integer, parameter :: NN=15, NSTACK=50
    real(dp):: a
    integer :: n,k,i,j,jstack,l,r
    integer, dimension(NSTACK) :: istack
    !-------
    n=size(arr)
    jstack=0
    l=1
    r=n
    do
      if (r-l < NN) then !Insertion sort when subarray small enough.
        do j=l+1,r
          a=arr(j)
          do i=j-1,l,-1
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a
        end do
        if (jstack==0) return
        r=istack(jstack) !________________Pop stack and begin a new round of partitioning
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+r)/2 !_______________________Choose median of left, center, and right elements as partitioning element a
        call SWAP_R(arr(k),arr(l+1)) !______Also rearrange so that a(l)<=a(l+1)<=a(r).
        call SWAP_RMASKED(arr(l),  arr(r),  arr(l)>arr(r))
        call SWAP_RMASKED(arr(l+1),arr(r),  arr(l+1)>arr(r))
        call SWAP_RMASKED(arr(l),  arr(l+1),arr(l)>arr(l+1))
        i=l+1
        j=r
        a=arr(l+1) !______________________Partitioning element
        do
          do !____________________________Scan up to find element >= a.
            i=i+1
            if(arr(i)>=a) exit
          end do
          do !____________________________Scan down to find element <= a
            j=j-1
            if(arr(j)<=a) exit
          end do
          if (j<i) exit !_______________Pointers crossed. Exit with partitioning complete.
          call SWAP_R(arr(i),arr(j))
        end do
        arr(l+1)=arr(j) !_________________Insert partitioning element.
        arr(j)=a
        jstack=jstack+2
        if (jstack > NSTACK) return !call nrerror('sort: NSTACK too small')
        if (r-i+1>=j-l) then
          istack(jstack)=r
          istack(jstack-1)=i
          r=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        end if
      end if
    end do
  end subroutine Sort

  !---

  subroutine CalcPermut(Arr,vPermut)
    !-----------------------------------------------------------------------------
    ! NR, Numerical Recipes for FORTRAN
    !-------------
    ! calc. permutation array vPermut that sorts Arr in increasing Arr(i)
    ! to apply permutation vPermut to array Arr: Arr=Arr(vPermut)
    !-----------------------------------------------------------------------------
    implicit none
    real(dp),dimension(:),intent(in) :: arr
    integer, dimension(:),intent(out):: vPermut
    integer, parameter:: NN=15, NSTACK=50
    real(dp) :: a
    integer :: n,k,i,j,indext,jstack,iLeft,iRite
    integer, dimension(NSTACK) :: istack
    !-------
    if(size(vPermut)/=size(arr)) return
    N=size(vPermut)
    do I=1,N
       vPermut(I)=I
    end do
    !vPermut=arth(1,1,n)
    JSTACK=0; iLeft=1; iRite=N
    do
       if (iRite-iLeft<NN) then
          do j=iLeft+1,iRite
             indext=vPermut(j)
             a=arr(indext)
             do i=j-1,iLeft,-1
                if (arr(vPermut(i)) <= a) exit
                vPermut(i+1)=vPermut(i)
             end do
             vPermut(i+1)=indext
          end do
          if (jstack==0) return
          iRite=istack(jstack)
          iLeft=istack(jstack-1)
          jstack=jstack-2
       else
          k=(iLeft+iRite)/2
          call SWAP_I(vPermut(k),vPermut(iLeft+1))
          call icomp_xchg(vPermut(iLeft),vPermut(iRite))
          call icomp_xchg(vPermut(iLeft+1),vPermut(iRite))
          call icomp_xchg(vPermut(iLeft),vPermut(iLeft+1))
          i=iLeft+1
          j=iRite
          indext=vPermut(iLeft+1)
          a=arr(indext)
          do
             do
                i=i+1
                if (arr(vPermut(i)) >= a) exit
             end do
             do
                j=j-1
                if (arr(vPermut(j)) <= a) exit
             end do
             if (j < i) exit
             call SWAP_I(vPermut(i),vPermut(j))
          end do
          vPermut(iLeft+1)=vPermut(j)
          vPermut(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) return !call nrerror('indexx: NSTACK too small')
          if (iRite-i+1 >= j-iLeft) then
             istack(jstack)=iRite
             istack(jstack-1)=i
             iRite=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=iLeft
             iLeft=i
          end if
       end if
    end do

  contains

    subroutine Icomp_Xchg(i,j) !in CalcPermut
      integer, intent(inout) :: i,j
      integer :: swp
      if (arr(j)<arr(i)) then
        swp=i; i=j; j=swp
      end if
    end subroutine Icomp_Xchg

  end subroutine CalcPermut

!--------------- FINITE DifFERENCE JACOBIAN ----------------------------
  subroutine Jacobian_Numeric(vFunc,X,vFuncX,tJac) !from "NR"
  !.forward-difference approximation to Jacobian
  !use M_Numeric_Tools
    real(dp),dimension(:),   intent(in)   :: vFuncX !vector(:) of function values at x
    real(dp),dimension(:),   intent(inout):: X !point(:) at which Jacobian is evaluated
    real(dp),dimension(:,:), intent(out)  :: tJac !(:,:) output Jacobian
    interface
       function vFunc(X)
         use M_Kinds
         implicit none
         real(dp),dimension(:),intent(in):: X
         real(dp),dimension(size(X))     :: vFunc
       end function vFunc
    endinterface
    !
    real(dp),parameter:: EPS=1.0E-6_dp
    real(dp),dimension(size(x)) :: Xsav,Xph,H
    integer:: J,N
    !
    N=   size(x)
    Xsav=X
    H=   EPS*abs(Xsav)
    where (H==Zero) H=Eps
    Xph=Xsav+H
    H=  Xph-Xsav !Trick to reduce finite precision error.
    do j=1,N
      X(J)= Xph(J)
      tJac(:,J)=(vFunc(x)-vFuncX(:))/H(J) !Forward difference formula.
      X(J)= Xsav(J)
    end do
  end subroutine Jacobian_Numeric

!-----------------------------------------------------------------------
  
  subroutine Jacobian_Numeric_Bis(vFunc,X,vFX,tJac) !from "NR"
    !.forward-difference approximation to Jacobian
    !use M_Numeric_Tools
    real(dp),dimension(:),   intent(in)   :: X      !point(:) at which Jacobian is evaluated
    real(dp),dimension(:),   intent(in)   :: vFX     !vector(:) of function values at x
    real(dp),dimension(:,:), intent(out)  :: tJac   !(:,:) output Jacobian
    interface
       function vFunc(X)
         use M_Kinds
         implicit none
         real(dp),dimension(:),intent(in):: X
         real(dp),dimension(size(X))     :: vFunc
       end function vFunc
    endinterface
    !
    real(dp),parameter:: EPS=1.0E-6_dp
    real(dp),dimension(size(x)) :: XplusH, H 
    integer:: J,N
    
    !--
    N = size(x)
    H = EPS*abs(X)

    do j=1,N
       if (H(j)==Zero) H(j)=Eps
    end do

    !--
    XplusH = X
    do j=1,N  
       XplusH(j) = X(J) + H(J)
       tJac(:,j) = ( vFunc(XplusH) - vFX )/ H(j) !Forward difference formula.
       XplusH(j) = X(j)
    end do

  end subroutine Jacobian_Numeric_Bis

subroutine linspace(xmin,xdelta,nval,vec)
  real(dp),intent(in) :: xmin,xdelta
  integer, intent(in) :: nval
  real(dp),intent(out):: vec(nval)
  !
  integer:: i
  !
  vec(1)= xmin
  do i=2,nval
    vec(i)= vec(1)+xdelta
  end do
end subroutine linspace

integer function bracket(m,xar,xi)
  integer, intent(in):: m       !the number of data values
  real(dp),intent(in):: xar(m)  !the data values, sorted
  real(dp),intent(in):: xi      !the query value
  !
  integer:: lo,mid,up
  !
  if(xi<xar(1) .or. xar(m)<xi) then
    bracket= -1
  else
    lo= 1    !lower
    up= mid  !upper
    do while(lo+1<up)
      mid= (lo+up)/2 !mid'point
      if (xi<xar(mid)) then ;  up= mid !mid'point becomes the new upper
      else                  ;  lo= mid !mid'point becomes the new lower
      end if
    end do
    bracket= lo
  end if
  !
end function bracket

subroutine linear_interp_2d(m,n,x1a,x2a,ya,x1,x2,  yo,ok)
  integer, intent(in) :: m,n
  real(dp),intent(in) :: x1a(m),x2a(n),ya(m,n)
  real(dp),intent(in) :: x1,x2
  real(dp),intent(out):: yo
  logical, intent(out):: ok
  !
  real(dp):: u,t
  integer :: j,k
  !
  ok= .true.
  j= bracket(m,x1a,x1)
  k= bracket(n,x2a,x2)
  !
  if(j<0 .or. k<0) then
    ok= .false.
    return
  end if
  !
  t= (x1 -x1a(j))/(x1a(j+1)-x1a(j))
  u= (x2 -x2a(j))/(x2a(k+1)-x2a(k))
  yo= (1-t)*(1-u)*ya(j,k)     &
  & +    t *(1-u)*ya(j+1,k)   &
  & +    t *u    *ya(j+1,k+1) &
  & + (1-t)*u    *ya(j,k+1)
  !
  return 
end subroutine linear_interp_2d

end module M_Numeric_Tools

