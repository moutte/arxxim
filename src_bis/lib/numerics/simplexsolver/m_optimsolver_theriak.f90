module M_Optimsolver_Theriak
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_
  implicit none
  !
  private
  !
  public:: Optimsolver_Theriak
  !
  real(dp), parameter:: XMin= 1.0D-9     !enforce vX(:)>Xmin
  real(dp), parameter:: DeltaMin= 1.0D-6 !min delta criteria in Theriak_Linesearch
  integer,  parameter:: m_max= 20        !max iter in Theriak_Linesearch
  !
  logical,parameter:: Debug= .true.  !! .false.  !! 
  integer,parameter:: IterMax= 100   !max iter in main loop
  !
  integer:: IterTot
  !
contains

subroutine Optimsolver_Theriak( & !
& ComputeMu,  & !
& ff,         & ! IN, log file index
& DeltaInit,  & ! IN
& TolX,       & ! IN
& vX,         & ! INOUT
& vMu,        & ! OUT
& G,          & ! OUT
& Converged,  & ! OUT
& Iter,nCallG)  ! OUT
  
  interface
    subroutine ComputeMu(vA,vB)
      use M_Kinds
      real(dp), intent(in) :: vA(:)
      real(dp), intent(out):: vB(:)
    end subroutine ComputeMu
  end interface
  
  integer, intent(in)   :: ff
  real(dp),intent(in)   :: DeltaInit
  real(dp),intent(in)   :: TolX
  real(dp),intent(inout):: vX(:)  ! composition
  real(dp),intent(out)  :: vMu(:) ! potentials
  real(dp),intent(out)  :: G      ! free energy
  logical, intent(out)  :: Converged
  integer, intent(out)  :: Iter,nCallG
  !---------------------------------------------------------------------
  real(dp):: vXlast1(size(vX)), vXlast2(size(vX))
  != former composition points
  real(dp):: vV(size(vX)),vVSteep(size(vX)),vXDif(size(vX))
  real(dp):: Gsav,DeltaG
  real(dp):: MuTilde,Slope
  real(dp):: Delta
  real(dp):: DeltaX
  integer :: N
  !---------------------------------------------------------------------
  if(iDebug>2) write(fTrc,'(/,A)') "< Optimsolver_Theriak"
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
  vXlast1(:)= vX(:)
  vXlast2(:)= vX(:)
  !
  if(FF>0) write(ff,'(/,A)') "========================"
  !
  Converged= .false.
  !
  nCallG= 0
  IterTot= 0
  Iter= 0
  call ComputeMu(vX,vMu)
  G= dot_product(vX,vMu)
  !
  do
  
    Iter= Iter+1
    
    !--- direction of steepest descent, normalized to |vVSteep|-1
    call ComputeMu(vX,vMu)
    MuTilde= SUM(vMu(:)) /real(N)
    vVSteep(:)= MuTilde -vMu(:)
    call Normalize(vVSteep)
    !---/
    
    if(Iter<3) then
      vV(:)= vVSteep(:)
    else
      vXDif(:)= vX(:) -vXlast2(:)
      Delta= Norm_L2(vXDif) *0.5D0
      Slope= dot_product(vXlast2(:),vMu(:)) - G
      !if(G > dot_product(vXlast2(:),vMu(:))) then
      if(Slope<Zero) then
        vV(:)= vVSteep(:)
      else
        vV(:)= vXDif(:) /Delta *0.5D0
      end if
      !print '(A,G12.3)',"Delta2=",Delta
    end if
    !
    vXlast2(:)= vXlast1(:)
    vXlast1(:)= vX(:)
    !
    Gsav= G
    call Theriak_Linesearch( &
    & ComputeMu,ff,Iter,N,Delta,vV, &
    & vX,vMu,G,nCallG)
    DeltaG= G-Gsav
    !
    if(Iter>2) then
      DeltaG= G-Gsav
      DeltaX= Norm_LInf(vX(:) -vXlast1(:))
      if(iDebug>2) write(fTrc,'(A,2G15.6)') "DeltaX_DeltaG=",DeltaX,DeltaG
      if(DeltaX < TolX) then
        if(iDebug>2) write(fTrc,'(A)') "==converged=YES=="
        Converged= .true.
        exit
      end if
    end if
    
    if(Iter>IterMax) then
      if(iDebug>2) write(fTrc,'(A,/)') "==converged=NO=="
      exit
    end if
  
  enddo
  !
  !print '(A,2F7.3)',"vXlast2",vXlast2(1),vXlast2(2)
  !print '(A,2F7.3)',"vXlast1",vXlast1(1),vXlast1(2)
  !print '(A,2F7.3)',"vX      ",vX(1),      vX(2)
  !
  if(iDebug>2) write(fTrc,'(A,/)') "</ Optimsolver_Theriak"
  
  return
end subroutine Optimsolver_Theriak

subroutine Theriak_Linesearch( &
& ComputeMu,ff,Iter,N,Delta,vV, &
& vX,vMu,G,nCallG)
  !---------------------------------------------------------------------
  interface
    subroutine ComputeMu(vA,vB)
      use M_Kinds
      real(dp), intent(in) :: vA(:)
      real(dp), intent(out):: vB(:)
    end subroutine ComputeMu
  end interface
  integer, intent(in):: ff
  integer, intent(in):: Iter
  integer, intent(in):: N
  real(dp),intent(in):: Delta
  real(dp),intent(in):: vV(:)
  real(dp),intent(inout):: vX(:)
  real(dp),intent(out)  :: vMu(:)
  real(dp),intent(inout):: G
  integer, intent(inout):: nCallG
  !---------------------------------------------------------------------
  integer:: m
  integer:: i
  real(dp):: G_Old,G_New,G_Last
  real(dp):: Delta_new
  character(len=80):: cForm
  !---------------------------------------------------------------------
  real(dp):: vXnew(N),vXlast(N)
  !real(dp):: vMu(N)
  !
  if(iDebug>2) write(cForm,'(a,i3,a)') '(A,A1,2(I3,A1),',N,'(G15.6,A1),G12.3)'
  !
  !---------------------------------------------- Additive linesearch --
  G_New= G
  do m=0,m_max
    !
    !--- increase step size
    !--- ( continue increase until G(vX +(m+1)*Delta)) > G(vX +m+*Delta) )
    Delta_new= (m+1)*Delta
    vXlast(:)= vXnew(:)
    G_Last= G_New
    !
    !--- compute new vX(:)
    vXnew(:)= vX(:) + vV(:) *Delta_new
    call Project(Xmin,vXnew)
    !
    G_Old= G_New
    call ComputeMu(vXnew,vMu)
    G_New= dot_product(vXnew,vMu)
    nCallG= nCallG +1
    !
    IterTot= IterTot+1
    if(FF>0) write(ff,cForm) "A",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
    !
    if(G_New > G_Old) exit
    !
  enddo
  !---------------------------------------------/ Additive linesearch --
  !
  !---------------------------------------------- Division linesearch --
  !--- ( case where G(vX+Delta)>G(vX) )
  if(m==0) then
    
    G_Old= G
    !
    do m= 0,m_max
      
      Delta_new= Delta /real(m+1)
      !
      vXnew(:)= vX(:) + vV(:) *Delta_new
      call Project(Xmin,vXnew)
      !
      vXlast(:)= vXnew(:)
      call ComputeMu(vXnew,vMu)
      G_New= dot_product(vXnew,vMu)
      nCallG= nCallG +1
      G_Last= G_New
      !
      IterTot= IterTot+1
      if(FF>0) write(ff,cForm) &
      & "B",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
      !
      if(        G_New < G_Old   &
      & .or. Delta_new < DeltaMin) exit
      
    enddo
    
  end if
  !---------------------------------------------/ Division linesearch --
  !
  vX(:)= vXlast(:)
  G= G_Last
  !
end subroutine Theriak_Linesearch

subroutine Normalize(vX)
  !--- vX - vX / Norm(vX)
  real(dp),intent(inout):: vX(:)
  real(dp):: NormvX
  !---
  NormvX= Norm_L2(vX)
  if(NormVX < 1.D-20 ) then
    vX= One
    NormvX= Norm_L2(vX)
  end if
  vX= vX /NormvX
  !
end subroutine Normalize
  
subroutine Project(Xminim,vX)
  !--
  !--- Project vX on {sum(vX) - 1, vX >-0}
  !--
  real(dp),intent(in)   :: Xminim
  real(dp),intent(inout):: vX(:)
  !
  real(dp):: SumV
  integer :: i
  !--
  !--- project on x >- Xminim
  do i=1,size(vX)
    if (vX(i)< Xminim) vX(i) = Xminim
  enddo
  !--- project on sum(xi) - 1 
  SumV= SUM(vX(:))
  vX(:)= vX(:) /SumV
  !
end subroutine Project

real(dp) function norm_L1(vX)
  !--- Compute Discrete L1 Norm
  real(dp),intent(in):: vX(:)
  !
  norm_L1 = SUM(ABS(vX(:)))
  !
end function norm_L1

real(dp) function norm_L2(vX)
  !-- -Compute Discrete L2 Norm
  real(dp),intent(in):: vX(:)
  !
  norm_L2= SQRT(dot_product(vX,vX) )
  !
end function norm_L2

real(dp) function norm_LInf(vX)
  !--- Compute Discrete LInfinity Norm
  real(dp),intent(in) :: vX(:)
  !
  norm_LInf= MAXVAL(ABS(vX))
  !
end function norm_LInf

end module M_Optimsolver_Theriak

! subroutine Optimsolver_Theriak_2( & !
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
! !-- and not as G- dot_product(vXnew,vMu)
! !-- -> Theriak_Linesearch_2 should be quicker
! !--
!   interface
!     real(dp) function Objective(vA)
!       use m_kinds
!       real(dp), intent(in) :: vA(:)
!     end function Objective
!     subroutine ComputeMu(X,vA,vB)
!       use M_Kinds
!       real(dp), intent(in) :: X
!       real(dp), intent(in) :: vA(:)
!       real(dp), intent(out):: vB(:)
!     end subroutine ComputeMu
!   end interface
!   integer, intent(in)   :: ff
!   real(dp),intent(in)   :: DeltaInit
!   real(dp),intent(in)   :: TolX
!   real(dp),intent(inout):: vX(:) ! composition
!   real(dp),intent(out)  :: G
!   logical, intent(out)  :: Converged
!   integer, intent(out)  :: Iter,nCallG
!   !
!   real(dp),allocatable:: vXlast1(:), vXlast2(:) ! former composition points
!   real(dp),allocatable:: vV(:),vVSteep(:)
!   real(dp),allocatable:: vXDif(:)
!   real(dp),allocatable:: vMu(:) !
!   real(dp):: Gsav,DeltaG
!   real(dp):: MuTilde,Slope
!   real(dp):: Delta
!   real(dp):: DeltaX
!   integer:: i,N
!   !
!   if(iDebug>2) write(fTrc,'(/,A)') "< Optimsolver_Theriak"
!   !
!   Delta= DeltaInit
!   !
!   N= size(vX)
!   !
!   allocate(vXlast1(N))
!   allocate(vXlast2(N))
!   allocate(vVSteep(N))
!   allocate(vV(N))
!   allocate(vXDif(N))
!   allocate(vMu(N))
!   !
!   vXlast1(:)= vX(:)
!   vXlast2(:)= vX(:)
!   !
!   if(FF>0) write(ff,'(/,A)') "========================"
!   !
!   Converged= .false.
!   !
!   nCallG= 0
!   IterTot= 0
!   Iter= 0
!   G= Objective(vX(:))
!   !
!   do
!     Iter= Iter+1
!     !
!     !--- direction of steepest descent, normalized to |vVSteep|-1
!     call ComputeMu(G,vX,vMu)
!     MuTilde= sum(vMu(:)) /real(N)
!     vVSteep(:)= MuTilde -vMu(:)
!     call Normalize(vVSteep)
!     !---/
!     !
!     if(Iter<3) then
!       vV(:)= vVSteep(:)
!     else
!       vXDif(:)= vX(:) -vXlast2(:)
!       Delta= Norm_L2(vXDif) *0.5D0
!       Slope= dot_product(vXlast2(:),vMu(:)) - G
!       !if(G > dot_product(vXlast2(:),vMu(:))) then
!       if(Slope<Zero) then
!         vV(:)= vVSteep(:)
!       else
!         vV(:)= vXDif(:) /Delta *0.5D0
!       end if
!       !print '(A,G12.3)',"Delta2=",Delta
!     end if
!     !
!     vXlast2(:)= vXlast1(:)
!     vXlast1(:)= vX(:)
!     !
!     Gsav= G
!     call Theriak_Linesearch_2(Objective,ff,Iter,N,Delta,vV,vX,G,nCallG)
!     DeltaG= G-Gsav
!     !
!     if(Iter>2) then
!       DeltaG= G-Gsav
!       DeltaX= norm_LInf(vX(:) -vXlast1(:))
!       if(iDebug>2) write(fTrc,'(A,2G15.6)') "DeltaX_DeltaG=",DeltaX,DeltaG
!       if(DeltaX < TolX) then
!         if(iDebug>2) write(fTrc,'(A)') "==converged=YES=="
!         Converged= .true.
!         exit
!       end if
!     end if
!     !
!     if(Iter>IterMax) then
!       if(iDebug>2) write(fTrc,'(A,/)') "==converged=NO=="
!       exit
!     end if
!   enddo
!   !
!   !print '(A,2F7.3)',"vXlast2",vXlast2(1),vXlast2(2)
!   !print '(A,2F7.3)',"vXlast1",vXlast1(1),vXlast1(2)
!   !print '(A,2F7.3)',"vX      ",vX(1),      vX(2)
!   !
!   deallocate(vXlast1)
!   deallocate(vXlast2)
!   deallocate(vVSteep)
!   deallocate(vV)
!   deallocate(vXDif)
!   deallocate(vMu)
!   !
!   if(iDebug>2) write(fTrc,'(A,/)') "</ Optimsolver_Theriak"
!   !
! end subroutine Optimsolver_Theriak_2
! 
! subroutine Theriak_Linesearch_2(Objective,ff,Iter,N,Delta,vV,vX,G,nCallG)
!   interface
!     real(dp) function Objective(vA)
!       use M_Kinds
!       real(dp),intent(in):: vA(:)
!     end function Objective
!   end interface
!   integer, intent(in):: ff
!   integer, intent(in):: Iter
!   integer, intent(in):: N
!   real(dp),intent(in):: Delta
!   real(dp),intent(in):: vV(:)
!   real(dp),intent(inout):: vX(:)
!   real(dp),intent(inout):: G
!   integer, intent(inout):: nCallG
!   !
!   integer:: m
!   integer:: i
!   real(dp):: G_Old,G_New,G_Last
!   real(dp):: Delta_new
!   character(len=80):: cForm
!   !
!   !real(dp),allocatable:: vXnew(:),vXlast(:)
!   real(dp):: vXnew(N),vXlast(N)
!   !
!   !allocate(vXnew(N))
!   !allocate(vXlast(N))
!   !
!   if(iDebug>2) write(cForm,'(a,i3,a)') '(A,A1,2(I3,A1),',N,'(G15.6,A1),G12.3)'
!   !
!   !---------------------------------------------- Additive linesearch --
!   G_New= G
!   do m=0,m_max
!     !
!     !--- increase step size
!     !--- ( continue increase until G(vX +(m+1)*Delta)) > G(vX +m+*Delta) )
!     Delta_new= (m+1)*Delta
!     vXlast(:)= vXnew(:)
!     G_Last= G_New
!     !
!     !--- compute new vX(:)
!     vXnew(:)= vX(:) + vV(:) *Delta_new
!     call Project(Xmin,vXnew)
!     !
!     G_Old= G_New
!     G_New= Objective(vXnew)  ;  nCallG= nCallG +1
!     !
!     IterTot= IterTot+1
!     if(FF>0) write(ff,cForm) "A",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
!     !
!     if(G_New > G_Old) exit
!     !
!   enddo
!   !---------------------------------------------/ Additive linesearch --
!   !
!   !---------------------------------------------- Division linesearch --
!   !--- ( case where G(vX+Delta)>G(vX) )
!   if(m==0) then
!     !
!     G_Old= G
!     !
!     do m= 0,m_max
!       !
!       Delta_new= Delta /real(m+1)
!       !
!       vXnew(:)= vX(:) + vV(:) *Delta_new
!       call Project(Xmin,vXnew)
!       !
!       vXlast(:)= vXnew(:)
!       G_New= Objective(vXnew)  ;  nCallG= nCallG +1
!       G_Last= G_New
!       !
!       IterTot= IterTot+1
!       if(FF>0) write(ff,cForm) "B",T_,Iter,T_,m,T_,(vXnew(i),T_,i=1,N),G_new
!       !
!       if(        G_New < G_Old   &
!       & .or. Delta_new < DeltaMin) exit
!     enddo
!   end if
!   !---------------------------------------------/ Division linesearch --
!   !
!   vX(:)= vXlast(:)
!   G= G_Last
!   !deallocate(vXnew)
!   !deallocate(vXlast)
!   !
! end subroutine Theriak_Linesearch_2
