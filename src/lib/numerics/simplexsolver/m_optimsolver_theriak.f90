MODULE M_Optimsolver_Theriak
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Optimsolver_Theriak
  !
  REAL(dp), PARAMETER:: XMin= 1.0D-9     !enforce vX(:)>Xmin
  REAL(dp), PARAMETER:: DeltaMin= 1.0D-6 !min delta criteria in Theriak_Linesearch
  INTEGER,  PARAMETER:: m_max= 20        !max iter in Theriak_Linesearch
  !
  LOGICAL,PARAMETER:: Debug= .TRUE.  !! .FALSE.  !! 
  INTEGER,PARAMETER:: IterMax= 100   !max iter in main loop
  !
  INTEGER:: IterTot
  !
CONTAINS

SUBROUTINE Optimsolver_Theriak( & !
& ComputeMu,  & !
& ff,         & ! IN, log file index
& DeltaInit,  & ! IN
& TolX,       & ! IN
& vX,         & ! INOUT
& vMu,        & ! OUT
& G,          & ! OUT
& Converged,  & ! OUT
& Iter,nCallG)  ! OUT
  
  INTERFACE
    SUBROUTINE ComputeMu(vA,vB)
      USE M_Kinds
      REAL(dp), INTENT(IN) :: vA(:)
      REAL(dp), INTENT(OUT):: vB(:)
    END SUBROUTINE ComputeMu
  END INTERFACE
  
  INTEGER, INTENT(IN)   :: ff
  REAL(dp),INTENT(IN)   :: DeltaInit
  REAL(dp),INTENT(IN)   :: TolX
  REAL(dp),INTENT(INOUT):: vX(:)  ! composition
  REAL(dp),INTENT(OUT)  :: vMu(:) ! potentials
  REAL(dp),INTENT(OUT)  :: G      ! free energy
  LOGICAL, INTENT(OUT)  :: Converged
  INTEGER, INTENT(OUT)  :: Iter,nCallG
  !---------------------------------------------------------------------
  REAL(dp),ALLOCATABLE:: vXlast1(:), vXlast2(:) ! former composition points
  REAL(dp),ALLOCATABLE:: vV(:),vVSteep(:)
  REAL(dp),ALLOCATABLE:: vXDif(:)
  !REAL(dp),ALLOCATABLE:: vMu(:) !
  REAL(dp):: Gsav,DeltaG
  REAL(dp):: MuTilde,Slope
  REAL(dp):: Delta
  REAL(dp):: DeltaX
  INTEGER :: N
  !---------------------------------------------------------------------
  IF(iDebug>2) WRITE(fTrc,'(/,A)') "< Optimsolver_Theriak"
  !
  !-- ALGORITHM THERIAK
  !--
  !-- start from arbitrary composition vX(:)
  !-- repeat
  !--   compute chemical potentials vMu(:) at vX(:)
  !--   compute steepest descent:
  !--     V(:)- sum(vMu(:))/N -Mu(:)
  !--     normalize: V(:)- V(:)/|V(:)|
  !--   adjust stepsize
  !--   compute new x(:)
  !-- until distance between successive x(:) is below tolerance
  !--
  !
  Delta= DeltaInit
  !
  N= size(vX)
  !
  ALLOCATE(vXlast1(N))
  ALLOCATE(vXlast2(N))
  ALLOCATE(vVSteep(N))
  ALLOCATE(vV(N))
  ALLOCATE(vXDif(N))
  !ALLOCATE(vMu(N))
  !
  vXlast1(:)= vX(:)
  vXlast2(:)= vX(:)
  !
  IF(FF>0) WRITE(ff,'(/,A)') "========================"
  !
  Converged= .FALSE.
  !
  nCallG= 0
  IterTot= 0
  Iter= 0
  CALL ComputeMu(vX,vMu)
  G= DOT_PRODUCT(vX,vMu)
  !
  DO
  
    Iter= Iter+1
    
    !--- direction of steepest descent, normalized to |vVSteep|-1
    CALL ComputeMu(vX,vMu)
    MuTilde= SUM(vMu(:)) /REAL(N)
    vVSteep(:)= MuTilde -vMu(:)
    CALL Normalize(vVSteep)
    !---/
    
    IF(Iter<3) THEN
      vV(:)= vVSteep(:)
    ELSE
      vXDif(:)= vX(:) -vXlast2(:)
      Delta= Norm_L2(vXDif) *0.5D0
      Slope= DOT_PRODUCT(vXlast2(:),vMu(:)) - G
      !IF(G > DOT_PRODUCT(vXlast2(:),vMu(:))) THEN
      IF(Slope<Zero) THEN
        vV(:)= vVSteep(:)
      ELSE
        vV(:)= vXDif(:) /Delta *0.5D0
      ENDIF
      !PRINT '(A,G12.3)',"Delta2=",Delta
    ENDIF
    !
    vXlast2(:)= vXlast1(:)
    vXlast1(:)= vX(:)
    !
    Gsav= G
    CALL Theriak_Linesearch( &
    & ComputeMu,ff,Iter,N,Delta,vV, &
    & vX,vMu,G,nCallG)
    DeltaG= G-Gsav
    !
    IF(Iter>2) THEN
      DeltaG= G-Gsav
      DeltaX= Norm_LInf(vX(:) -vXlast1(:))
      IF(iDebug>2) WRITE(fTrc,'(A,2G15.6)') "DeltaX_DeltaG=",DeltaX,DeltaG
      IF(DeltaX < TolX) THEN
        IF(iDebug>2) WRITE(fTrc,'(A)') "==converged=YES=="
        Converged= .TRUE.
        EXIT
      ENDIF
    ENDIF
    
    IF(Iter>IterMax) THEN
      IF(iDebug>2) WRITE(fTrc,'(A,/)') "==converged=NO=="
      EXIT
    ENDIF
  
  ENDDO
  !
  !PRINT '(A,2F7.3)',"vXlast2",vXlast2(1),vXlast2(2)
  !PRINT '(A,2F7.3)',"vXlast1",vXlast1(1),vXlast1(2)
  !PRINT '(A,2F7.3)',"vX      ",vX(1),      vX(2)
  !
  DEALLOCATE(vXlast1)
  DEALLOCATE(vXlast2)
  DEALLOCATE(vVSteep)
  DEALLOCATE(vV)
  DEALLOCATE(vXDif)
  !DEALLOCATE(vMu)
  !
  IF(iDebug>2) WRITE(fTrc,'(A,/)') "</ Optimsolver_Theriak"
  
  RETURN
END SUBROUTINE Optimsolver_Theriak

SUBROUTINE Theriak_Linesearch( &
& ComputeMu,ff,Iter,N,Delta,vV, &
& vX,vMu,G,nCallG)
  !---------------------------------------------------------------------
  INTERFACE
    SUBROUTINE ComputeMu(vA,vB)
      USE M_Kinds
      REAL(dp), INTENT(IN) :: vA(:)
      REAL(dp), INTENT(OUT):: vB(:)
    END SUBROUTINE ComputeMu
  END INTERFACE
  INTEGER, INTENT(IN):: ff
  INTEGER, INTENT(IN):: Iter
  INTEGER, INTENT(IN):: N
  REAL(dp),INTENT(IN):: Delta
  REAL(dp),INTENT(IN):: vV(:)
  REAL(dp),INTENT(INOUT):: vX(:)
  REAL(dp),INTENT(OUT)  :: vMu(:)
  REAL(dp),INTENT(INOUT):: G
  INTEGER, INTENT(INOUT):: nCallG
  !---------------------------------------------------------------------
  INTEGER:: m
  INTEGER:: i
  REAL(dp):: G_Old,G_New,G_Last
  REAL(dp):: Delta_new
  CHARACTER(LEN=80):: cForm
  !---------------------------------------------------------------------
  REAL(dp):: vXnew(N),vXlast(N)
  !REAL(dp):: vMu(N)
  !
  IF(iDebug>2) WRITE(cForm,'(a,i3,a)') '(A,A1,2(I3,A1),',N,'(G15.6,A1),G12.3)'
  !
  !---------------------------------------------- Additive linesearch --
  G_New= G
  DO m=0,m_max
    !
    !--- increase step size
    !--- ( continue increase until G(vX +(m+1)*Delta)) > G(vX +m+*Delta) )
    Delta_new= (m+1)*Delta
    vXlast(:)= vXnew(:)
    G_Last= G_New
    !
    !--- compute new vX(:)
    vXnew(:)= vX(:) + vV(:) *Delta_new
    CALL Project(Xmin,vXnew)
    !
    G_Old= G_New
    CALL ComputeMu(vXnew,vMu)
    G_New= DOT_PRODUCT(vXnew,vMu)
    nCallG= nCallG +1
    !
    IterTot= IterTot+1
    IF(FF>0) WRITE(ff,cForm) "A",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
    !
    IF(G_New > G_Old) EXIT
    !
  ENDDO
  !---------------------------------------------/ Additive linesearch --
  !
  !---------------------------------------------- Division linesearch --
  !--- ( case where G(vX+Delta)>G(vX) )
  IF(m==0) THEN
    
    G_Old= G
    !
    DO m= 0,m_max
      
      Delta_new= Delta /REAL(m+1)
      !
      vXnew(:)= vX(:) + vV(:) *Delta_new
      CALL Project(Xmin,vXnew)
      !
      vXlast(:)= vXnew(:)
      CALL ComputeMu(vXnew,vMu)
      G_New= DOT_PRODUCT(vXnew,vMu)
      nCallG= nCallG +1
      G_Last= G_New
      !
      IterTot= IterTot+1
      IF(FF>0) WRITE(ff,cForm) &
      & "B",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
      !
      IF(        G_New < G_Old   &
      & .OR. Delta_new < DeltaMin) EXIT
      
    ENDDO
    
  ENDIF
  !---------------------------------------------/ Division linesearch --
  !
  vX(:)= vXlast(:)
  G= G_Last
  !
END SUBROUTINE Theriak_Linesearch

SUBROUTINE Normalize(vX)
  !--- vX - vX / Norm(vX)
  REAL(dp),INTENT(INOUT):: vX(:)
  REAL(dp):: NormvX
  !---
  NormvX= Norm_L2(vX)
  IF(NormVX < 1.D-20 ) THEN
    vX= One
    NormvX= Norm_L2(vX)
  end if
  vX= vX /NormvX
  !
END SUBROUTINE Normalize
  
SUBROUTINE Project(Xminim,vX)
  !--
  !--- Project vX on {sum(vX) - 1, vX >-0}
  !--
  REAL(dp),INTENT(IN)   :: Xminim
  REAL(dp),INTENT(INOUT):: vX(:)
  !
  REAL(dp):: SumV
  INTEGER :: i
  !--
  !--- project on x >- Xminim
  DO i=1,SIZE(vX)
    IF (vX(i)< Xminim) vX(i) = Xminim
  ENDDO
  !--- project on sum(xi) - 1 
  SumV= SUM(vX(:))
  vX(:)= vX(:) /SumV
  !
END SUBROUTINE Project

REAL(dp) FUNCTION norm_L1(vX)
  !--- Compute Discrete L1 Norm
  REAL(dp),INTENT(IN):: vX(:)
  !
  norm_L1 = SUM(ABS(vX(:)))
  !
END FUNCTION norm_L1

REAL(dp) FUNCTION norm_L2(vX)
  !-- -Compute Discrete L2 Norm
  REAL(dp),INTENT(IN):: vX(:)
  !
  norm_L2= SQRT(DOT_PRODUCT(vX,vX) )
  !
END FUNCTION norm_L2

REAL(dp) FUNCTION norm_LInf(vX)
  !--- Compute Discrete LInfinity Norm
  REAL(dp),INTENT(IN) :: vX(:)
  !
  norm_LInf= MAXVAL(ABS(vX))
  !
END FUNCTION norm_LInf

END MODULE M_Optimsolver_Theriak

! SUBROUTINE Optimsolver_Theriak_2( & !
! & Objective,  & !
! & ComputeMu,  & !
! & ff,         & ! IN, log file index
! & DeltaInit,  & ! IN
! & TolX,       & ! IN
! & vX,         & ! INOUT
! & G,          & ! OUT
! & Converged,  & ! OUT
! & Iter,nCallG)  ! OUT
! !--
! !-- this version use routine Objective to compute Gibbs directly
! !-- and not as G- DOT_PRODUCT(vXnew,vMu)
! !-- -> Theriak_Linesearch_2 should be quicker
! !--
!   INTERFACE
!     REAL(dp) FUNCTION Objective(vA)
!       USE m_kinds
!       REAL(dp), INTENT(IN) :: vA(:)
!     END FUNCTION Objective
!     SUBROUTINE ComputeMu(X,vA,vB)
!       USE M_Kinds
!       REAL(dp), INTENT(IN) :: X
!       REAL(dp), INTENT(IN) :: vA(:)
!       REAL(dp), INTENT(OUT):: vB(:)
!     END SUBROUTINE ComputeMu
!   END INTERFACE
!   INTEGER, INTENT(IN)   :: ff
!   REAL(dp),INTENT(IN)   :: DeltaInit
!   REAL(dp),INTENT(IN)   :: TolX
!   REAL(dp),INTENT(INOUT):: vX(:) ! composition
!   REAL(dp),INTENT(OUT)  :: G
!   LOGICAL, INTENT(OUT)  :: Converged
!   INTEGER, INTENT(OUT)  :: Iter,nCallG
!   !
!   REAL(dp),ALLOCATABLE:: vXlast1(:), vXlast2(:) ! former composition points
!   REAL(dp),ALLOCATABLE:: vV(:),vVSteep(:)
!   REAL(dp),ALLOCATABLE:: vXDif(:)
!   REAL(dp),ALLOCATABLE:: vMu(:) !
!   REAL(dp):: Gsav,DeltaG
!   REAL(dp):: MuTilde,Slope
!   REAL(dp):: Delta
!   REAL(dp):: DeltaX
!   INTEGER:: i,N
!   !
!   IF(iDebug>2) WRITE(fTrc,'(/,A)') "< Optimsolver_Theriak"
!   !
!   Delta= DeltaInit
!   !
!   N= size(vX)
!   !
!   ALLOCATE(vXlast1(N))
!   ALLOCATE(vXlast2(N))
!   ALLOCATE(vVSteep(N))
!   ALLOCATE(vV(N))
!   ALLOCATE(vXDif(N))
!   ALLOCATE(vMu(N))
!   !
!   vXlast1(:)= vX(:)
!   vXlast2(:)= vX(:)
!   !
!   IF(FF>0) WRITE(ff,'(/,A)') "========================"
!   !
!   Converged= .FALSE.
!   !
!   nCallG= 0
!   IterTot= 0
!   Iter= 0
!   G= Objective(vX(:))
!   !
!   DO
!     Iter= Iter+1
!     !
!     !--- direction of steepest descent, normalized to |vVSteep|-1
!     CALL ComputeMu(G,vX,vMu)
!     MuTilde= sum(vMu(:)) /REAL(N)
!     vVSteep(:)= MuTilde -vMu(:)
!     CALL Normalize(vVSteep)
!     !---/
!     !
!     IF(Iter<3) THEN
!       vV(:)= vVSteep(:)
!     ELSE
!       vXDif(:)= vX(:) -vXlast2(:)
!       Delta= Norm_L2(vXDif) *0.5D0
!       Slope= DOT_PRODUCT(vXlast2(:),vMu(:)) - G
!       !IF(G > DOT_PRODUCT(vXlast2(:),vMu(:))) THEN
!       IF(Slope<Zero) THEN
!         vV(:)= vVSteep(:)
!       ELSE
!         vV(:)= vXDif(:) /Delta *0.5D0
!       ENDIF
!       !PRINT '(A,G12.3)',"Delta2=",Delta
!     ENDIF
!     !
!     vXlast2(:)= vXlast1(:)
!     vXlast1(:)= vX(:)
!     !
!     Gsav= G
!     CALL Theriak_Linesearch_2(Objective,ff,Iter,N,Delta,vV,vX,G,nCallG)
!     DeltaG= G-Gsav
!     !
!     IF(Iter>2) THEN
!       DeltaG= G-Gsav
!       DeltaX= norm_LInf(vX(:) -vXlast1(:))
!       IF(iDebug>2) WRITE(fTrc,'(A,2G15.6)') "DeltaX_DeltaG=",DeltaX,DeltaG
!       IF(DeltaX < TolX) THEN
!         IF(iDebug>2) WRITE(fTrc,'(A)') "==converged=YES=="
!         Converged= .TRUE.
!         EXIT
!       ENDIF
!     ENDIF
!     !
!     IF(Iter>IterMax) THEN
!       IF(iDebug>2) WRITE(fTrc,'(A,/)') "==converged=NO=="
!       EXIT
!     ENDIF
!   ENDDO
!   !
!   !PRINT '(A,2F7.3)',"vXlast2",vXlast2(1),vXlast2(2)
!   !PRINT '(A,2F7.3)',"vXlast1",vXlast1(1),vXlast1(2)
!   !PRINT '(A,2F7.3)',"vX      ",vX(1),      vX(2)
!   !
!   DEALLOCATE(vXlast1)
!   DEALLOCATE(vXlast2)
!   DEALLOCATE(vVSteep)
!   DEALLOCATE(vV)
!   DEALLOCATE(vXDif)
!   DEALLOCATE(vMu)
!   !
!   IF(iDebug>2) WRITE(fTrc,'(A,/)') "</ Optimsolver_Theriak"
!   !
! END SUBROUTINE Optimsolver_Theriak_2
! 
! SUBROUTINE Theriak_Linesearch_2(Objective,ff,Iter,N,Delta,vV,vX,G,nCallG)
!   INTERFACE
!     REAL(dp) FUNCTION Objective(vA)
!       USE M_Kinds
!       REAL(dp),INTENT(IN):: vA(:)
!     END FUNCTION Objective
!   END INTERFACE
!   INTEGER, INTENT(IN):: ff
!   INTEGER, INTENT(IN):: Iter
!   INTEGER, INTENT(IN):: N
!   REAL(dp),INTENT(IN):: Delta
!   REAL(dp),INTENT(IN):: vV(:)
!   REAL(dp),INTENT(INOUT):: vX(:)
!   REAL(dp),INTENT(INOUT):: G
!   INTEGER, INTENT(INOUT):: nCallG
!   !
!   INTEGER:: m
!   INTEGER:: i
!   REAL(dp):: G_Old,G_New,G_Last
!   REAL(dp):: Delta_new
!   CHARACTER(LEN=80):: cForm
!   !
!   !real(dp),allocatable:: vXnew(:),vXlast(:)
!   REAL(dp):: vXnew(N),vXlast(N)
!   !
!   !allocate(vXnew(N))
!   !allocate(vXlast(N))
!   !
!   IF(iDebug>2) WRITE(cForm,'(a,i3,a)') '(A,A1,2(I3,A1),',N,'(G15.6,A1),G12.3)'
!   !
!   !---------------------------------------------- Additive linesearch --
!   G_New= G
!   DO m=0,m_max
!     !
!     !--- increase step size
!     !--- ( continue increase until G(vX +(m+1)*Delta)) > G(vX +m+*Delta) )
!     Delta_new= (m+1)*Delta
!     vXlast(:)= vXnew(:)
!     G_Last= G_New
!     !
!     !--- compute new vX(:)
!     vXnew(:)= vX(:) + vV(:) *Delta_new
!     CALL Project(Xmin,vXnew)
!     !
!     G_Old= G_New
!     G_New= Objective(vXnew)  ;  nCallG= nCallG +1
!     !
!     IterTot= IterTot+1
!     IF(FF>0) WRITE(ff,cForm) "A",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
!     !
!     IF(G_New > G_Old) EXIT
!     !
!   ENDDO
!   !---------------------------------------------/ Additive linesearch --
!   !
!   !---------------------------------------------- Division linesearch --
!   !--- ( case where G(vX+Delta)>G(vX) )
!   IF(m==0) THEN
!     !
!     G_Old= G
!     !
!     DO m= 0,m_max
!       !
!       Delta_new= Delta /REAL(m+1)
!       !
!       vXnew(:)= vX(:) + vV(:) *Delta_new
!       CALL Project(Xmin,vXnew)
!       !
!       vXlast(:)= vXnew(:)
!       G_New= Objective(vXnew)  ;  nCallG= nCallG +1
!       G_Last= G_New
!       !
!       IterTot= IterTot+1
!       IF(FF>0) WRITE(ff,cForm) "B",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
!       !
!       IF(        G_New < G_Old   &
!       & .OR. Delta_new < DeltaMin) EXIT
!     ENDDO
!   ENDIF
!   !---------------------------------------------/ Division linesearch --
!   !
!   vX(:)= vXlast(:)
!   G= G_Last
!   !deallocate(vXnew)
!   !deallocate(vXlast)
!   !
! END SUBROUTINE Theriak_Linesearch_2
