MODULE M_Numeric_Tools 

  !---------------------------------------------------------
  ! "NR" Tools from Numerical Recipes"
  !---------------------------------------------------------

  USE M_Kinds
  USE M_Trace, ONLY:Stop_
  IMPLICIT NONE
 
  PRIVATE
  
  !---
  PUBLIC:: Interpol
  !---
  PUBLIC:: vAbs
  !---
  PUBLIC:: Swap_RV
  PUBLIC:: Swap_I
  !---
  PUBLIC:: OuterProd_R
  PUBLIC:: OuterAnd
  !---
  PUBLIC:: iMaxLoc_R
  PUBLIC:: iMinLoc_R
  PUBLIC:: iFirstLoc
  !---
  PUBLIC:: Unit_Matrix
  PUBLIC:: Sort
  PUBLIC:: CalcPermut
  !---
  PUBLIC:: Jacobian_Numeric
  PUBLIC:: Jacobian_Numeric_Bis
  !---

  INTEGER,PUBLIC:: fNewt_I= 0 !counter used for DebNewt
  INTEGER,PUBLIC:: fNewtF=  0 !"trace" file for Newton, values of the unknown vector at each step, opened if DebNewt
  INTEGER,PUBLIC:: fNewtR=  0 !"trace" file for Newton, values of the residue at each step, opened if DebNewt 
  INTEGER,PUBLIC:: fNewtJ=  0 !
  !
  !INTEGER,PARAMETER::NPAR_ARTH=16,NPAR2_ARTH=8 !for ARTH_I

CONTAINS


 !--------------- LINEAR TABLE INTERPOLATION  -------------------------------------------


  REAL(dp) FUNCTION Interpol(T,iota,N,vT,vP)
    !----------------------------------------------------------------------
    ! provisional, basic linear interpolation of P from a table (vT,vP)
    ! assumes the vT-series is strictly increasing from 1 to N
    !----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp),              INTENT(IN)   :: T,iota
    INTEGER,               INTENT(IN)   :: N
    REAL(dp), DIMENSION(:),INTENT(IN)   :: vT,vP
    !
    REAL(dp):: P,K1,K2
    INTEGER :: I
    !
    IF(T<vT(1) .OR. T>vT(N)) CALL Stop_("T outside T-series")
    !
    I=1; DO WHILE(T>vT(I)); I=I+1; ENDDO !-> find upper bound
    !
    IF(ABS(T-vT(I))<=Iota) THEN
       P= vP(I)
    ELSE
       K1= (T-vT(I)  ) /(vT(I-1) - vT(I)  )
       K2= (T-vT(I-1)) /(vT(I)   - vT(I-1))
       P=  K1 *vP(I-1) + K2 *vP(I)
    ENDIF
    Interpol= P
    RETURN
  ENDFUNCTION Interpol


  !--------------- VECTOR UTILS  -------------------------------------------


  REAL(dp) FUNCTION vAbs(V) 
    !----------------------
    ! Norm of Vector V
    !----------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:),INTENT(IN) :: V
    vAbs=SQRT(DOT_PRODUCT(V,V))
  ENDFUNCTION vAbs

  !---

  SUBROUTINE Swap_RV(A,B) 
    !-----------------------
    ! Swap RealVector
    !-----------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:),INTENT(INOUT):: A,B
    REAL(dp),DIMENSION(SIZE(A))        :: DUM
    DUM=A; A=B; B=DUM
  ENDSUBROUTINE Swap_RV

  !---

  FUNCTION iMaxLoc_R(ARR) 
    !-------------------------------
    ! Localize maximum value
    !-------------------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:), INTENT(IN) :: ARR
    INTEGER               :: iMaxLoc_R
    INTEGER, DIMENSION(1) :: IMAX
    IMAX=MAXLOC(ARR(:))
    iMaxLoc_R=IMAX(1)
  END FUNCTION iMaxLoc_R

  !---

  FUNCTION iMinLoc_R(ARR) 
    !-------------------------------
    ! Localize minimum value
    !-------------------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:), INTENT(IN) :: ARR
    INTEGER               :: iMinLoc_R
    INTEGER, DIMENSION(1) :: IMin
    IMin=MINLOC(ARR(:))
    iMinLoc_R=IMin(1)
  END FUNCTION iMinLoc_R

  !---

  FUNCTION iFirstLoc(MASK) 
    !-------------------------------
    ! Localize first occurence
    !-------------------------------
    IMPLICIT NONE
    INTEGER                        :: iFirstLoc
    LOGICAL,DIMENSION(:),INTENT(IN):: MASK
    !
    INTEGER,DIMENSION(1):: LOC
    !
    LOC=MAXLOC(MERGE(1,0,MASK))
    iFirstLoc=LOC(1)
    IF (.NOT. MASK(IFIRSTLOC)) IFIRSTLOC=SIZE(MASK)+1
  ENDFUNCTION IFIRSTLOC

  !---

  FUNCTION OuterProd_R(A,B) 
    !-------------------------------
    ! Outer product of vectors
    !-------------------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:), INTENT(IN) :: A,B
    REAL(dp),DIMENSION(SIZE(A),SIZE(B)) :: OuterProd_R
    !
    OuterProd_R= SPREAD(A,DIM=2,NCOPIES=SIZE(B)) &
         &          * SPREAD(B,DIM=1,NCOPIES=SIZE(A))
  END FUNCTION OuterProd_R

  !---

  FUNCTION OuterAND(A,B) 
    !-------------------------------
    ! Outer operator and 
    !-------------------------------
    IMPLICIT NONE
    LOGICAL, DIMENSION(:), INTENT(IN)   :: A,B
    LOGICAL, DIMENSION(SIZE(A),SIZE(B)) :: OuterAND
    OuterAND= SPREAD(A,DIM=2,NCOPIES=SIZE(B)) &
         &   .AND. SPREAD(B,DIM=1,NCOPIES=SIZE(A))
  END FUNCTION OuterAND

  !---

  SUBROUTINE Unit_Matrix(MAT) 
    !------------------------------------------------------------
    ! Square real matrix Identity, dimension=MIN(nRows,nColumns)
    !------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:,:),INTENT(OUT) :: MAT
    INTEGER:: I,N
    N=MIN(SIZE(MAT,1),SIZE(MAT,2))
    MAT(:,:)=Zero
    DO I=1,N
       MAT(I,I)=One
    ENDDO
  ENDSUBROUTINE UNIT_MATRIX


  !--------------- SORTING ---------------------------------------------------

  SUBROUTINE SWAP_R(A,B) 
    !-----------------------
    ! Swap Real
    !-----------------------
    IMPLICIT NONE
    REAL(dp),INTENT(INOUT) :: A,B
    REAL(dp) :: DUM
    DUM=A; A=B; B=DUM
  ENDSUBROUTINE SWAP_R

  !---

  SUBROUTINE SWAP_I(A,B)
    !-----------------------
    ! Swap Integer
    !-----------------------
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: A,B
    INTEGER:: DUM
    DUM=A; A=B; B=DUM
  ENDSUBROUTINE SWAP_I

  !---

  SUBROUTINE SWAP_RMASKED(A,B,MASK)
    !-----------------------
    ! Swap Real with Mask
    !-----------------------
    IMPLICIT NONE
    REAL(dp),INTENT(INOUT) :: A,B
    LOGICAL,  INTENT(IN)   :: MASK
    REAL(dp):: SWP
    IF (MASK) THEN
       SWP=A; A=B; B=SWP
    ENDIF
  ENDSUBROUTINE SWAP_RMASKED

  !---

  SUBROUTINE Sort(arr) 
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
    IMPLICIT NONE
    REAL(dp),DIMENSION(:), INTENT(INOUT) :: arr
    INTEGER, PARAMETER :: NN=15, NSTACK=50
    REAL(dp):: a
    INTEGER :: n,k,i,j,jstack,l,r
    INTEGER, DIMENSION(NSTACK) :: istack
    !-------
    n=SIZE(arr)
    jstack=0
    l=1
    r=n
    DO
       IF (r-l < NN) THEN !Insertion sort when subarray small enough.
          DO j=l+1,r
             a=arr(j)
             DO i=j-1,l,-1
                IF (arr(i) <= a) EXIT
                arr(i+1)=arr(i)
             ENDDO
             arr(i+1)=a
          ENDDO
          IF (jstack==0) RETURN
          r=istack(jstack) !________________Pop stack and begin a new round of partitioning
          l=istack(jstack-1)
          jstack=jstack-2
       ELSE
          k=(l+r)/2 !_______________________Choose median of left, center, and right elements as partitioning element a
          CALL SWAP_R(arr(k),arr(l+1)) !______Also rearrange so that a(l)<=a(l+1)<=a(r).
          CALL SWAP_RMASKED(arr(l),  arr(r),  arr(l)>arr(r))
          CALL SWAP_RMASKED(arr(l+1),arr(r),  arr(l+1)>arr(r))
          CALL SWAP_RMASKED(arr(l),  arr(l+1),arr(l)>arr(l+1))
          i=l+1
          j=r
          a=arr(l+1) !______________________Partitioning element
          DO
             DO !____________________________Scan up to find element >= a.
                i=i+1
                IF(arr(i)>=a) EXIT
             ENDDO
             DO !____________________________Scan down to find element <= a
                j=j-1
                IF(arr(j)<=a) EXIT
             ENDDO
             IF (j<i) EXIT !_______________Pointers crossed. Exit with partitioning complete.
             CALL SWAP_R(arr(i),arr(j))
          ENDDO
          arr(l+1)=arr(j) !_________________Insert partitioning element.
          arr(j)=a
          jstack=jstack+2
          IF (jstack > NSTACK) RETURN !CALL nrerror('sort: NSTACK too small')
          IF (r-i+1>=j-l) THEN
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          ELSE
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          ENDIF
       ENDIF
    ENDDO
  ENDSUBROUTINE Sort

  !---

  SUBROUTINE CalcPermut(Arr,vPermut)
    !-----------------------------------------------------------------------------
    ! NR, Numerical Recipes for FORTRAN
    !-------------
    ! calc. permutation array vPermut that sorts Arr in increasing Arr(i)
    ! to apply permutation vPermut to array Arr: Arr=Arr(vPermut)
    !-----------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp),DIMENSION(:),INTENT(IN) :: arr
    INTEGER, DIMENSION(:),INTENT(OUT):: vPermut
    INTEGER, PARAMETER:: NN=15, NSTACK=50
    REAL(dp) :: a
    INTEGER :: n,k,i,j,indext,jstack,iLeft,iRite
    INTEGER, DIMENSION(NSTACK) :: istack
    !-------
    IF(SIZE(vPermut)/=SIZE(arr)) RETURN
    N=SIZE(vPermut)
    DO I=1,N
       vPermut(I)=I
    ENDDO
    !vPermut=arth(1,1,n)
    JSTACK=0; iLeft=1; iRite=N
    DO
       IF (iRite-iLeft<NN) THEN
          DO j=iLeft+1,iRite
             indext=vPermut(j)
             a=arr(indext)
             DO i=j-1,iLeft,-1
                IF (arr(vPermut(i)) <= a) EXIT
                vPermut(i+1)=vPermut(i)
             ENDDO
             vPermut(i+1)=indext
          ENDDO
          IF (jstack==0) RETURN
          iRite=istack(jstack)
          iLeft=istack(jstack-1)
          jstack=jstack-2
       ELSE
          k=(iLeft+iRite)/2
          CALL SWAP_I(vPermut(k),vPermut(iLeft+1))
          CALL icomp_xchg(vPermut(iLeft),vPermut(iRite))
          CALL icomp_xchg(vPermut(iLeft+1),vPermut(iRite))
          CALL icomp_xchg(vPermut(iLeft),vPermut(iLeft+1))
          i=iLeft+1
          j=iRite
          indext=vPermut(iLeft+1)
          a=arr(indext)
          DO
             DO
                i=i+1
                IF (arr(vPermut(i)) >= a) EXIT
             ENDDO
             DO
                j=j-1
                IF (arr(vPermut(j)) <= a) EXIT
             ENDDO
             IF (j < i) EXIT
             CALL SWAP_I(vPermut(i),vPermut(j))
          ENDDO
          vPermut(iLeft+1)=vPermut(j)
          vPermut(j)=indext
          jstack=jstack+2
          IF (jstack > NSTACK) RETURN !CALL nrerror('indexx: NSTACK too small')
          IF (iRite-i+1 >= j-iLeft) THEN
             istack(jstack)=iRite
             istack(jstack-1)=i
             iRite=j-1
          ELSE
             istack(jstack)=j-1
             istack(jstack-1)=iLeft
             iLeft=i
          ENDIF
       ENDIF
    ENDDO

  CONTAINS

    SUBROUTINE Icomp_Xchg(i,j) !in CalcPermut
      INTEGER, INTENT(INOUT) :: i,j
      INTEGER :: swp
      IF (arr(j)<arr(i)) THEN
         swp=i; i=j; j=swp
      ENDIF
    ENDSUBROUTINE Icomp_Xchg

  ENDSUBROUTINE CalcPermut


  !--------------- FINITE DIFFERENCE JACOBIAN -------------------------------------------
  

  SUBROUTINE Jacobian_Numeric(vFunc,X,vFuncX,tJac) !from "NR"
    !.forward-difference approximation to Jacobian
    !USE M_Numeric_Tools
    REAL(dp),DIMENSION(:),   INTENT(IN)   :: vFuncX !vector(:) of function values at x
    REAL(dp),DIMENSION(:),   INTENT(INOUT):: X !point(:) at which Jacobian is evaluated
    REAL(dp),DIMENSION(:,:), INTENT(OUT)  :: tJac !(:,:) output Jacobian
    INTERFACE
       FUNCTION vFunc(X)
         USE M_Kinds
         IMPLICIT NONE
         REAL(dp),DIMENSION(:),INTENT(IN):: X
         REAL(dp),DIMENSION(SIZE(X))     :: vFunc
       ENDFUNCTION vFunc
    ENDINTERFACE
    !
    REAL(dp),PARAMETER:: EPS=1.0E-6_dp
    REAL(dp),DIMENSION(SIZE(x)) :: Xsav,Xph,H
    INTEGER:: J,N
    !
    N=   SIZE(x)
    Xsav=X
    H=   EPS*ABS(Xsav)
    WHERE (H==Zero) H=Eps
    Xph=Xsav+H
    H=  Xph-Xsav !Trick to reduce finite precision error.
    DO j=1,N
      X(J)= Xph(J)
      tJac(:,J)=(vFunc(x)-vFuncX(:))/H(J) !Forward difference formula.
      X(J)= Xsav(J)
    ENDDO
  ENDSUBROUTINE Jacobian_Numeric

  !---

  SUBROUTINE Jacobian_Numeric_Bis(vFunc,X,vFX,tJac) !from "NR"
    !.forward-difference approximation to Jacobian
    !USE M_Numeric_Tools
    REAL(dp),DIMENSION(:),   INTENT(IN)   :: X      !point(:) at which Jacobian is evaluated
    REAL(dp),DIMENSION(:),   INTENT(IN)   :: vFX     !vector(:) of function values at x
    REAL(dp),DIMENSION(:,:), INTENT(OUT)  :: tJac   !(:,:) output Jacobian
    INTERFACE
       FUNCTION vFunc(X)
         USE M_Kinds
         IMPLICIT NONE
         REAL(dp),DIMENSION(:),INTENT(IN):: X
         REAL(dp),DIMENSION(SIZE(X))     :: vFunc
       ENDFUNCTION vFunc
    ENDINTERFACE
    !
    REAL(dp),PARAMETER:: EPS=1.0E-6_dp
    REAL(dp),DIMENSION(SIZE(x)) :: XplusH, H 
    INTEGER:: J,N
    
    !--
    N = SIZE(x)
    H = EPS*ABS(X)

    DO j=1,N
       if (H(j)==Zero) H(j)=Eps
    end DO

    !--
    XplusH = X
    DO j=1,N  
       XplusH(j) = X(J) + H(J)
       tJac(:,j) = ( vFunc(XplusH) - vFX )/ H(j) !Forward difference formula.
       XplusH(j) = X(j)
    ENDDO

  ENDSUBROUTINE Jacobian_Numeric_Bis

ENDMODULE M_Numeric_Tools

