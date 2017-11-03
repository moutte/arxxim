module M_Basis_Ortho
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_, Stop_
  implicit none
  !
  private
  !
  public:: Basis_Ortho_Compute
  public:: Basis_Ortho_Compute_Debug

  !! public:: OrthoBasis
  !! public:: Basis_Ortho_Alloc
  !! public:: Basis_Ortho_Clear
  !
  !! real(dp),dimension(:,:),allocatable:: tOrthoAlfSp

  !! public :: tOrthoAlfSp

contains

  subroutine Basis_Ortho_Compute(vIndexBuffer,tAlfSp,tXOrthoBasis)
    integer, intent(in) :: vIndexBuffer(:)
    real(dp),intent(in) :: tAlfSp(:,:)
    real(dp),intent(out):: tXOrthoBasis(:,:)
    !
    real(dp), allocatable :: tVBasis(:,:)
    real(dp), allocatable :: tEBasis(:,:) !, tXOrthoBasis(:,:)
    real(dp), allocatable :: tOrthoBasis(:,:), tNormBasis(:)
    integer :: nE, nV, nCi, nCp
    integer :: iE, iV, iCp
    logical :: Ok
    ! real(dp) :: coef
    !---

    if (idebug>1) write(fTrc,'(/,A)') "< Basis_Ortho_Compute"

    !// Allocate Working Arrays

    nE =  size(vIndexBuffer)
    nV =  count(vIndexBuffer(:)>0)
    nCp=  nE
    nCi=  nE - nV

    allocate( tVBasis(nV, nE) )          ! MOBILE or BUFFER COMPONENTS
    allocate( tEBasis(nE, nE) )          ! PRIMARY COMPONENTS
    allocate( tOrthoBasis(nE, nE) )
    allocate( tNormBasis(nE) )
    !allocate( tXOrthoBasis(nE, nE) )

    iV= 0
    do iCp=1,nCp
      if(vIndexBuffer(iCp)>0) then
        iV= iV +1
        tVBasis(iV,:)= tAlfSp(:,vIndexBuffer(iCp))
      end if
    end do

    do iE = 1, nE
      tEBasis(iE,:)  = 0.D0
      tEBasis(iE,iE) = 1.D0
    end do

    !// Compute Basis

    call OrthoBasis( &
    & 0,tVBasis, tEBasis, nV, nE, &
    & tOrthoBasis, tNormBasis, Ok)
    if ( .not.( Ok ) ) call Stop_("Error OrthoBasis")

    do iCp = 1, nCi
      iV = nV + iCp
      tXOrthoBasis ( iCp, 1:nE ) = tOrthoBasis(iV,1:nE)
    end do
    do iCp = nCi+1, nCp
      iV = iCp - nCi
      tXOrthoBasis ( iCp, 1:nE ) = tOrthoBasis(iV,1:nE)
    end do

    if(iDebug>2) call Matrix_Sho(tOrthoBasis)

    !// Compute Formula Matrices

    !tOrthoAlfSp  = matmul( tXOrthoBasis , tAlfSp  )

    !// Deallocate Working Arrays

    deallocate ( tVBasis )
    deallocate ( tEBasis )
    deallocate ( tOrthoBasis )
    deallocate ( tNormBasis )

    if (idebug>1) write(fTrc,'(A,/)') "</ Basis_Ortho_Compute"

  end subroutine Basis_Ortho_Compute

  subroutine Matrix_Sho(T)
    use M_IoTools,    only: GetUnit
    use M_Files,      only: DirOut

    real(dp),dimension(:,:),intent(in):: T
    !
    integer:: nS,nC
    integer:: F
    integer:: I,J
    !
    call GetUnit(F)
    open(F,file=trim(DirOut)//"_othobasis.log")
    !call Files_Index_Write(fHtm,...
    !
    nC= size(T,1)
    nS= size(T,2)
    !
    do I=1,nC
      do J=1,nS
        write(F,'(F7.3,A1)',advance="no") T(I,J), T_
      end do
      write(F,*)
    end do
    !
    close(F)
    !
  end subroutine Matrix_Sho

  subroutine Basis_Ortho_Compute_Debug(vIndexBuffer,tAlfSp,tXOrthoBasis)
    use M_IoTools,    only: GetUnit
    use M_Files,      only: DirOut
    !
    integer,          intent(in) :: vIndexBuffer(:)
    real(dp),         intent(in) :: tAlfSp(:,:)
    real(dp),         intent(out):: tXOrthoBasis(:,:)
    !
    real(dp), allocatable :: tVBasis(:, :)
    real(dp), allocatable :: tEBasis(:, :)
    real(dp), allocatable :: tOrthoBasis(:,:), tNormBasis(:)
    integer :: nE, nV, nCi, nCp
    integer :: iE, iV, J, iCp
    integer :: F
    logical :: Ok
    real(dp) :: coef
    !---

    call GetUnit(F)
    open(F,file=trim(DirOut)//"_orthobasis.log")

    if (idebug>1) write(F,*) "Basis_Ortho_Compute"
    !
    !// Allocate Working Arrays
    nCp= size(vIndexBuffer)
    nE = nCp
    nV = count(vIndexBuffer(:)>0)

    if( nV==0) return !-------------------------------------------return

    nCi= nE - nV

    allocate(tVBasis(nV,nE))          ! MOBILE or BUFFER COMPONENTS
    allocate(tEBasis(nE,nE))          ! PRIMARY COMPONENTS
    allocate(tOrthoBasis(nE,nE))
    allocate(tNormBasis(nE))

    !!do iV = 1, nV
    !!  J = nCi + iV
    !!  tVBasis(iV,:)= tAlfSp(:,vCpn(J)%iSpc)
    !!end do

    iV= 0
    do iE=1,nE
      if(vIndexBuffer(iE)>0) then
        iV= iV +1
        tVBasis(iV,:)= tAlfSp(:,vIndexBuffer(iE))
      end if
    end do
    do iE = 1, nE
      tEBasis(iE,:)  = 0.D0
      tEBasis(iE,iE) = 1.D0
    end do

    !// Compute Basis

    if (idebug>1) then
      write(F,'(A)') "================================================"
      write(F,'(A)') 'Compute Orthogonal Basis'
      write(F,'(A)') "================================================"
    end if

    call OrthoBasis( &
    & F,tVBasis, tEBasis, nV, nE, &
    & tOrthoBasis, tNormBasis, Ok)
    if ( .not.( Ok ) ) call Stop_("Error OrthoBasis")

    !// Show Ortho Basis Results

    if (idebug>1) then
      write(F,'(A)') "================================================"
      write(F,'(A)') 'Orthogonal Formulas System = tOrthoBasis'
      write(F,'(A)') "================================================"
    end if

    if (idebug>1) then
      do iV = 1, nE
        if (iV<=nV) write(F,'(A,I3,A)',advance="no") 'Buffer Basis       ', iV, ' = '
        if (iV>nV)  write(F,'(A,I3,A)',advance="no") 'Orthogonal Formula ', iV-nV, ' = '
        do iE = 1, nE
          coef = tOrthoBasis(iV, iE)
          !! if ( abs(coef) > 0.00001) write(F,'(A4,A,F8.3,A)',advance="no") &
            !! & trim( vEle(iCpnEle(iE))%NamEl),'(',coef, ') '
          if ( abs(coef) > 0.00001) write(F,'(I3,A,F8.3,A)',advance="no") &
            & iE,'(',coef, ') '
        end do

        write(F,*) ' '
        if ( iV == nV ) write(F,*) "----------"
      end do
    end if

    !// Reorder Equations For Arxim Dynamic Problem

    if (idebug>1) then
      write(F,'(A)') "================================================"
      write(F,'(A)') 'Mole Balance Equations Re-Ordered for Arxim = tXOrthoBasis'
      write(F,'(A)') "================================================"
    end if
    
    do iCp= 1, nCi
      iV= nV + iCp
      tXOrthoBasis(iCp,1:nE) = tOrthoBasis(iV,1:nE)
    end do
    do iCp= nCi+1, nCp
      iV= iCp - nCi
      tXOrthoBasis(iCp,1:nE) = tOrthoBasis(iV,1:nE)
    end do

    if (idebug>1) then
      do iV= 1, nE
        !! if (iV<=nCi)  write(F,'(A,I3,A,A,A)',advance="no") 'Balance Equation Ci', iV ,      '  [', vEle(iCpnEle(iV)) % NamEl, '] = '
        !! if (iV>nCi)   write(F,'(A,I3,A,A,A)',advance="no") 'Balance Buffer   Bx', iV - nCi, '  [', vEle(iCpnEle(iV)) % NamEl, '] = '
        if (iV<=nCi) write(F,'(A,I3,A)',advance="no") &
        & 'Balance Equation Ci', iV ,      ' = '
        if (iV>nCi)  write(F,'(A,I3,A)',advance="no") &
        & 'Balance Buffer   Bx', iV - nCi, ' = '
        !
        do iE= 1, nE
          coef = tXOrthoBasis(iV, iE)
          !! if ( abs(coef) > 0.00001) write(F,'(A4,A,F8.3,A)',advance="no") &
          !! & trim( vEle(iCpnEle(iE))%NamEl),'(',coef, ') '
          if ( abs(coef) > 0.00001) write(F,'(I3,A,F8.3,A)',advance="no") &
          & iE,'(',coef, ') '
        end do
        !
        write(F,*) ' '
        if (iV==nCi) write(F,*) "----------"
      end do
    end if

    !// Compute Formula Matrices

    if (idebug>1) then
      write(F,'(A)') "================================================"
      write(F,'(A)') 'Compute MODIFIED Formula Matrices'
      write(F,'(A)') "================================================"
      write(F,*) 'tOrthoAlfSp  = matmul( tXOrthoBasis , tAlfSp  )'
    end if

    !tOrthoAlfSp  = matmul( tXOrthoBasis , tAlfSp  )

    !// Deallocate Working Arrays

    deallocate (tVBasis)
    deallocate (tEBasis)
    deallocate (tOrthoBasis)
    deallocate (tNormBasis)
    !deallocate ( tXOrthoBasis )

    if (idebug>1) then
      write(F,'(A)') "================================================"
      write(F,'(A)') "Basis_Ortho_Basis Compute Done !"
    end if

    close(F)

  end subroutine Basis_Ortho_Compute_Debug

  !----

  subroutine OrthoBasis(FF,tVBasis, tEBasis, nV, nE, tOrthoBasis, tNormBasis, Ok)
    !-------------------------------------------
    ! Purpose : Orthogonalize Basis of Subspace
    !           Using Gramm Schmitt Method
    !-------------------------------------------

    implicit none

    integer, intent(in) :: FF                         ! trace file
    integer, intent(in) :: nV                         ! Dimension of SubSpace V
    integer, intent(in) :: nE                         ! Dimension of Space V
    real(dp), intent(in)  :: tVBasis(nV, nE)          ! Basis of E
    real(dp), intent(in)  :: tEBasis(nE, nE)          ! Basis of E
    real(dp), intent(out) :: tOrthoBasis(nE, nE)      ! Ortho Basis of E = [V|V']
    real(dp), intent(out) :: tNormBasis(nE)           ! Norm of Ortho Basis of E = [V|V']
    logical, intent(out) :: Ok
    !-----
    real(dp) :: CheckZero = 1.d-6
    real(dp) :: vj(nE), vk(nE)
    real(dp) :: normvk, pjk

    integer :: iV, iE, k, j

    if (FF>0) write(FF,'(A,2I3,A)') "OrthoBasis : nE,nV =[", nE,  nV, ']'

    Ok = .false.
    tOrthoBasis = 0.D0
    tNormBasis  = 1.D0

    !------------
    ! Compute a basis of tVBasis
    k = 1


    do iV = 1, nV

      !write(21,*) "k=", k

      vk = tVBasis(iV,:)

      do j = 1, k-1
        vj     = tOrthoBasis(j,:)
        pjk    = dot_product (vk,vj)
        vk     = vk - pjk*vj
      end do

      normvk = sqrt(dot_product(vk,vk))
      !write(21,*) "iV =", iV, " Normvk = ", normvK

      if (normvk < CheckZero ) then
        stop "VBasis Problem"
      else
        vk = vk/normvk
        if (FF>0) write(FF,'(A,I3,A,20F8.3)') "[",k,"]", vk(1:nE)
        if ( k > nE ) stop "OrthoBasis Max ESize Problem"
        tOrthoBasis(k,:) = vk
        tNormBasis(k)    = normvk
        k = k + 1
      end if

    end do

    !------------
    ! Complete by vectors of tEBasis
    k = nV + 1

    do iE =1, nE

      !write(21,*) "k=", k


      vk = tEBasis(iE,:)

      do j = 1, k-1
        vj     = tOrthoBasis(j,:)
        pjk    = dot_product (vk,vj)
        vk     = vk - pjk*vj
      end do

      normvk = sqrt(dot_product(vk,vk))
      !write(21,*) "iE =", iE, " Normvk = ", normvK

      if (normvk < CheckZero) then
        !write(21,*) "Vector", iE, "is Zero Projected"
      else
        vk = vk/normvk
        if (FF>0) write(FF,'(A,I3,A,20F8.3)') "[",k,"]", vk(1:nE)
        if ( k > nE ) stop "OrthoBasis Max ESize Problem"
        tOrthoBasis(k,:) = vk
        tNormBasis(k)    = normvk
        k = k + 1
      end if

    end do

    if ( k <= nE ) stop "Total ESize Problem"

    Ok = .true.

  end subroutine OrthoBasis

  !! subroutine Basis_Ortho_Alloc(tAlfSp)
  !! real(dp),intent(in):: tAlfSp(:,:)
  !! !---
  !! !if (idebug>1) write(21,*) "Basis_Ortho_Alloc"

  !! allocate(tOrthoAlfSp(size(tAlfSp,1), size(tAlfSp,2)))

  !! if (idebug>1) then
  !! write(21,'(A)') "=========================================================="
  !! write(21,*) "Shape of tOrthoAlfSp  = ",  SHAPE(tOrthoAlfSp)
  !! write(21,'(A)') "=========================================================="
  !! end if

  !! end subroutine Basis_Ortho_Alloc

  !! !---

  !! subroutine Basis_Ortho_Clear

  !! deallocate(tOrthoAlfSp)

  !! end subroutine Basis_Ortho_Clear

  !---

end module M_Basis_Ortho
