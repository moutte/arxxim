MODULE M_Basis_Ortho
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_, Stop_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Basis_Ortho_Compute
  PUBLIC:: Basis_Ortho_Compute_Debug

  !~ PUBLIC:: OrthoBasis
  !~ PUBLIC:: Basis_Ortho_Alloc
  !~ PUBLIC:: Basis_Ortho_Clear
  !
  !~ REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tOrthoAlfSp

  !~ PUBLIC :: tOrthoAlfSp

CONTAINS

  SUBROUTINE Basis_Ortho_Compute(vIndexBuffer,tAlfSp,tXOrthoBasis)
    INTEGER, INTENT(IN) :: vIndexBuffer(:)
    REAL(dp),INTENT(IN) :: tAlfSp(:,:)
    REAL(dp),INTENT(OUT):: tXOrthoBasis(:,:)
    !
    REAL(dp), ALLOCATABLE :: tVBasis(:, :)
    REAL(dp), ALLOCATABLE :: tEBasis(:, :) !, tXOrthoBasis(:,:)
    REAL(dp), ALLOCATABLE :: tOrthoBasis(:,:), tNormBasis(:)
    INTEGER :: nE, nV, nCi, nCp
    INTEGER :: iE, iV, iCp
    LOGICAL :: Ok
    ! REAL(dp) :: coef
    !---

    IF (iDebug>0) write(fTrc,'(/,A)') "< Basis_Ortho_Compute"

    !// Allocate Working Arrays

    nE =  SIZE(vIndexBuffer)
    nV =  COUNT(vIndexBuffer(:)>0)
    nCp=  nE
    nCi=  nE - nV

    ALLOCATE( tVBasis(nV, nE) )          ! MOBILE or BUFFER COMPONENTS
    ALLOCATE( tEBasis(nE, nE) )          ! PRIMARY COMPONENTS
    ALLOCATE( tOrthoBasis(nE, nE) )
    ALLOCATE( tNormBasis(nE) )
    !ALLOCATE( tXOrthoBasis(nE, nE) )

    iV= 0
    DO iCp=1,nCp
      IF(vIndexBuffer(iCp)>0) THEN
        iV= iV +1
        tVBasis(iV,:)= tAlfSp(:,vIndexBuffer(iCp))
      ENDIF
    END DO

    DO iE = 1, nE
      tEBasis(iE,:)  = 0.D0
      tEBasis(iE,iE) = 1.D0
    END DO

    !// Compute Basis

    CALL OrthoBasis( &
    & 0,tVBasis, tEBasis, nV, nE, &
    & tOrthoBasis, tNormBasis, Ok)
    IF ( .not.( Ok ) ) CALL Stop_("Error OrthoBasis")

    DO iCp = 1, nCi
      iV = nV + iCp
      tXOrthoBasis ( iCp, 1:nE ) = tOrthoBasis(iV,1:nE)
    END DO
    DO iCp = nCi+1, nCp
      iV = iCp - nCi
      tXOrthoBasis ( iCp, 1:nE ) = tOrthoBasis(iV,1:nE)
    END DO

    IF(iDebug>2) CALL Matrix_Sho(tOrthoBasis)

    !// Compute Formula Matrices

    !tOrthoAlfSp  = MATMUL( tXOrthoBasis , tAlfSp  )

    !// Deallocate Working Arrays

    DEALLOCATE ( tVBasis )
    DEALLOCATE ( tEBasis )
    DEALLOCATE ( tOrthoBasis )
    DEALLOCATE ( tNormBasis )
    !DEALLOCATE ( tXOrthoBasis )

    IF (iDebug>0) write(fTrc,'(A,/)') "</ Basis_Ortho_Compute"

  END SUBROUTINE Basis_Ortho_Compute

  SUBROUTINE Matrix_Sho(T)
    USE M_IoTools,    ONLY: GetUnit
    USE M_Files,      ONLY: DirOut

    REAL(dp),DIMENSION(:,:),INTENT(IN):: T
    !
    INTEGER:: nS,nC
    INTEGER:: F
    INTEGER:: I,J
    !
    CALL GetUnit(F)
    OPEN(F,FILE=TRIM(DirOut)//"_othobasis.log")
    !CALL Files_Index_Write(fHtm,...
    !
    nC= SIZE(T,1)
    nS= SIZE(T,2)
    !
    DO I=1,nC
      DO J=1,nS
        WRITE(F,'(F7.3,A1)',ADVANCE='NO') T(I,J), T_
      ENDDO
      WRITE(F,*)
    ENDDO
    !
    CLOSE(F)
    !
  END SUBROUTINE Matrix_Sho

  SUBROUTINE Basis_Ortho_Compute_Debug(vIndexBuffer,tAlfSp,tXOrthoBasis)
    USE M_IoTools,    ONLY: GetUnit
    USE M_Files,      ONLY: DirOut
    USE M_Global_Vars,ONLY: vEle
    USE M_System_Vars,ONLY: vCpn
    !
    INTEGER, INTENT(IN) :: vIndexBuffer(:)
    REAL(dp),INTENT(IN) :: tAlfSp(:,:)
    REAL(dp),INTENT(OUT):: tXOrthoBasis(:,:)
    !
    REAL(dp), ALLOCATABLE :: tVBasis(:, :)
    REAL(dp), ALLOCATABLE :: tEBasis(:, :)
    REAL(dp), ALLOCATABLE :: tOrthoBasis(:,:), tNormBasis(:)
    INTEGER :: nE, nV, nCi, nCp
    INTEGER :: iE, iV, J, iCp
    INTEGER :: F
    LOGICAL :: Ok
    REAL(dp) :: coef
    !---

    CALL GetUnit(F)
    OPEN(F,FILE=TRIM(DirOut)//"_orthobasis.log")

    IF (iDebug>0) write(F,*) "Basis_Ortho_Compute"
    !
    !// Allocate Working Arrays

    nCp= SIZE(vCpn)
    nCi= COUNT(vCpn(:)%Statut=="INERT")

    nE = nCp
    nV = COUNT(vIndexBuffer(:)>0)

    ALLOCATE( tVBasis(nV, nE) )          ! MOBILE or BUFFER COMPONENTS
    ALLOCATE( tEBasis(nE, nE) )          ! PRIMARY COMPONENTS
    ALLOCATE( tOrthoBasis(nE, nE) )
    ALLOCATE( tNormBasis(nE) )
    !ALLOCATE( tXOrthoBasis(nE, nE) )

    DO iV = 1, nV
      J = nCi + iV
      tVBasis(iV,:)= tAlfSp(:,vCpn(J)%iSpc)
    END DO

    !~ iV= 0
    !~ DO iE=1,nE
    !~ IF(vIndexBuffer(iE)>0) THEN
    !~ iV= iV +1
    !~ tVBasis(iV,:)= tAlfSp(:,vIndexBuffer(iE))
    !~ ENDIF
    !~ END DO

    DO iE = 1, nE
      tEBasis(iE,:)  = 0.D0
      tEBasis(iE,iE) = 1.D0
    END DO

    !// Compute Basis

    IF (iDebug>0) THEN
      WRITE(F,'(A)') "================================================"
      WRITE(F,'(A)') 'Compute Orthogonal Basis'
      WRITE(F,'(A)') "================================================"
    END IF

    CALL OrthoBasis( &
    & F,tVBasis, tEBasis, nV, nE, &
    & tOrthoBasis, tNormBasis, Ok)
    if ( .not.( Ok ) ) call Stop_("Error OrthoBasis")

    !// Show Ortho Basis Results

    IF (iDebug>0) THEN
      WRITE(F,'(A)') "=========================================================="
      WRITE(F,'(A)') 'Orthogonal Formulas System = tOrthoBasis'
      WRITE(F,'(A)') "=========================================================="
    END IF

    IF (iDebug>0) THEN
      DO iV = 1, nE
        IF (iV<=nV) WRITE(F,'(A,I3,A)',ADVANCE='NO') 'Buffer Basis       ', iV, ' = '
        IF (iV>nV)  WRITE(F,'(A,I3,A)',ADVANCE='NO') 'Orthogonal Formula ', iV-nV, ' = '
        DO iE = 1, nE
          coef = tOrthoBasis(iV, iE)
          !~ IF ( abs(coef) > 0.00001) WRITE(F,'(A4,A,F8.3,A)',ADVANCE='NO') &
            !~ & TRIM( vEle(iCpnEle(iE))%NamEl),'(',coef, ') '
          IF ( abs(coef) > 0.00001) WRITE(F,'(I3,A,F8.3,A)',ADVANCE='NO') &
            & iE,'(',coef, ') '
        END DO

        WRITE(F,*) ' '
        IF ( iV == nV ) WRITE(F,*) "----------"
      END DO
    END IF

    !// Reorder Equations For Arxim Dynamic Problem

    IF (iDebug>0) THEN
      WRITE(F,'(A)') "=========================================================="
      WRITE(F,'(A)') 'Mole Balance Equations Re-Ordered for Arxim = tXOrthoBasis'
      WRITE(F,'(A)') "=========================================================="
    END IF

    DO iCp = 1, nCi
      iV = nV + iCp
      tXOrthoBasis ( iCp, 1:nE ) = tOrthoBasis(iV,1:nE)
    END DO
    DO iCp = nCi+1, nCp
      iV = iCp - nCi
      tXOrthoBasis ( iCp, 1:nE ) = tOrthoBasis(iV,1:nE)
    END DO

    IF (iDebug>0) THEN
      DO iV = 1, nE
        !~ IF (iV<=nCi)  WRITE(F,'(A,I3,A,A,A)',ADVANCE='NO') 'Balance Equation Ci', iV ,      '  [', vEle(iCpnEle(iV)) % NamEl, '] = '
        !~ IF (iV>nCi)   WRITE(F,'(A,I3,A,A,A)',ADVANCE='NO') 'Balance Buffer   Bx', iV - nCi, '  [', vEle(iCpnEle(iV)) % NamEl, '] = '
        IF (iV<=nCi) WRITE(F,'(A,I3,A)',ADVANCE='NO') &
        & 'Balance Equation Ci', iV ,      ' = '
        IF (iV>nCi)  WRITE(F,'(A,I3,A)',ADVANCE='NO') &
        & 'Balance Buffer   Bx', iV - nCi, ' = '

        DO iE = 1, nE
          coef = tXOrthoBasis(iV, iE)
          !~ IF ( abs(coef) > 0.00001) write(F,'(A4,A,F8.3,A)',ADVANCE='NO') &
            !~ & TRIM( vEle(iCpnEle(iE))%NamEl),'(',coef, ') '
          IF ( abs(coef) > 0.00001) write(F,'(I3,A,F8.3,A)',ADVANCE='NO') &
            & iE,'(',coef, ') '
        END DO

        WRITE(F,*) ' '
        IF ( iV == nCi ) WRITE(F,*) "----------"
      END DO
    END IF

    !// Compute Formula Matrices

    IF (iDebug>0) THEN
      WRITE(F,'(A)') "====================================================="
      WRITE(F,'(A)') 'Compute Modified Formula Matrices'
      WRITE(F,'(A)') "====================================================="
      WRITE(F,*) 'tOrthoAlfSp  = matmul( tXOrthoBasis , tAlfSp  )'
    END IF

    !tOrthoAlfSp  = MATMUL( tXOrthoBasis , tAlfSp  )

    !// Deallocate Working Arrays

    DEALLOCATE ( tVBasis )
    DEALLOCATE ( tEBasis )
    DEALLOCATE ( tOrthoBasis )
    DEALLOCATE ( tNormBasis )
    !DEALLOCATE ( tXOrthoBasis )

    IF (iDebug>0) THEN
      WRITE(F,'(A)') "====================================================="
      WRITE(F,'(A)') "Basis_Ortho_Basis Compute Done !"
    END IF

    CLOSE(F)

  END SUBROUTINE Basis_Ortho_Compute_Debug

  !----

  SUBROUTINE OrthoBasis(FF,tVBasis, tEBasis, nV, nE, tOrthoBasis, tNormBasis, Ok)
    !-------------------------------------------
    ! Purpose : Orthogonalize Basis of Subspace
    !           Using Gramm Schmitt Method
    !-------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: FF                         ! trace file
    INTEGER, INTENT(in) :: nV                         ! Dimension of SubSpace V
    INTEGER, INTENT(in) :: nE                         ! Dimension of Space V
    REAL(dp), INTENT(in)  :: tVBasis(nV, nE)          ! Basis of E
    REAL(dp), INTENT(in)  :: tEBasis(nE, nE)          ! Basis of E
    REAL(dp), INTENT(out) :: tOrthoBasis(nE, nE)      ! Ortho Basis of E = [V|V']
    REAL(dp), INTENT(out) :: tNormBasis(nE)           ! Norm of Ortho Basis of E = [V|V']
    LOGICAL, INTENT(out) :: Ok
    !-----
    REAL(dp) :: CheckZero = 1.d-6
    REAL(dp) :: vj(nE), vk(nE)
    REAL(dp) :: normvk, pjk

    INTEGER :: iV, iE, k, j

    IF (FF>0) WRITE(FF,'(A,2I3,A)') "OrthoBasis : nE,nV =[", nE,  nV, ']'

    Ok = .FALSE.
    tOrthoBasis = 0.D0
    tNormBasis  = 1.D0

    !------------
    ! Compute a basis of tVBasis
    k = 1


    DO iV = 1, nV

      !write(21,*) "k=", k

      vk = tVBasis(iV,:)

      DO j = 1, k-1
        vj     = tOrthoBasis(j,:)
        pjk    = DOT_PRODUCT (vk,vj)
        vk     = vk - pjk*vj
      END DO

      normvk = SQRT(DOT_PRODUCT(vk,vk))
      !write(21,*) "iV =", iV, " Normvk = ", normvK

      IF (normvk < CheckZero ) THEN
        STOP "VBasis Problem"
      ELSE
        vk = vk/normvk
        IF (FF>0) WRITE(FF,'(A,I3,A,20F8.3)') "[",k,"]", vk(1:nE)
        IF ( k > nE ) STOP "OrthoBasis Max ESize Problem"
        tOrthoBasis(k,:) = vk
        tNormBasis(k)    = normvk
        k = k + 1
      END IF

    END DO

    !------------
    ! Complete by vectors of tEBasis
    k = nV + 1

    DO iE =1, nE

      !write(21,*) "k=", k


      vk = tEBasis(iE,:)

      DO j = 1, k-1
        vj     = tOrthoBasis(j,:)
        pjk    = DOT_PRODUCT (vk,vj)
        vk     = vk - pjk*vj
      END DO

      normvk = SQRT(DOT_PRODUCT(vk,vk))
      !write(21,*) "iE =", iE, " Normvk = ", normvK

      IF (normvk < CheckZero) THEN
        !write(21,*) "Vector", iE, "is Zero Projected"
      ELSE
        vk = vk/normvk
        IF (FF>0) WRITE(FF,'(A,I3,A,20F8.3)') "[",k,"]", vk(1:nE)
        IF ( k > nE ) STOP "OrthoBasis Max ESize Problem"
        tOrthoBasis(k,:) = vk
        tNormBasis(k)    = normvk
        k = k + 1
      END IF

    END DO

    IF ( k <= nE ) STOP "Total ESize Problem"

    Ok = .TRUE.

  END SUBROUTINE OrthoBasis

  !~ SUBROUTINE Basis_Ortho_Alloc(tAlfSp)
  !~ REAL(dp),INTENT(IN):: tAlfSp(:,:)
  !~ !---
  !~ !IF (iDebug>0) write(21,*) "Basis_Ortho_Alloc"

  !~ ALLOCATE(tOrthoAlfSp(SIZE(tAlfSp,1), SIZE(tAlfSp,2)))

  !~ IF (iDebug>0) THEN
  !~ WRITE(21,'(A)') "=========================================================="
  !~ WRITE(21,*) "Shape of tOrthoAlfSp  = ",  SHAPE(tOrthoAlfSp)
  !~ WRITE(21,'(A)') "=========================================================="
  !~ END IF

  !~ END SUBROUTINE Basis_Ortho_Alloc

  !~ !---

  !~ SUBROUTINE Basis_Ortho_Clear

  !~ DEALLOCATE(tOrthoAlfSp)

  !~ END SUBROUTINE Basis_Ortho_Clear

  !---

END MODULE M_Basis_Ortho
