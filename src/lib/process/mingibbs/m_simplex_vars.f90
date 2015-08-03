MODULE M_Simplex_Vars
!--
!-- vars for simplex computations
!--
  USE M_Kinds
  USE M_T_MixModel,ONLY: MaxPole
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Simplex_Vars_Alloc
  PUBLIC:: Simplex_Vars_Clean
  PUBLIC:: tSimplex,IZROV,IPOSV
  !
  REAL(dp),ALLOCATABLE:: tSimplex(:,:)
  INTEGER, ALLOCATABLE:: IZROV(:),IPOSV(:)

CONTAINS

SUBROUTINE Simplex_Vars_Alloc(M,N)
  INTEGER,INTENT(IN)::M,N
  !
  CALL Simplex_Vars_Clean
  !
  ALLOCATE(tSimplex(0:M+1,0:N))  ;  tSimplex=Zero
  ALLOCATE(IPOSV(1:M))           ;  IPOSV= 0
  ALLOCATE(IZROV(1:N))           ;  IZROV= 0
  !
ENDSUBROUTINE Simplex_Vars_Alloc

SUBROUTINE Simplex_Vars_Clean
  IF(ALLOCATED(tSimplex))  DEALLOCATE(tSimplex)
  IF(ALLOCATED(IZROV))     DEALLOCATE(IZROV)
  IF(ALLOCATED(IPOSV))     DEALLOCATE(IPOSV)
ENDSUBROUTINE Simplex_Vars_Clean

END MODULE M_Simplex_Vars
