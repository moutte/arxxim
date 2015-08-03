MODULE M_T_TPCond
!--
!-- T,P conditions
!--

  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: T_TPCond
  !
  !! PUBLIC:: TPcond_IndexTP
  !! PUBLIC:: TPcond_IndexT
  !
  TYPE:: T_TPCond !container for T,P Properties
    CHARACTER(15)::Name
    REAL(dp):: TdgC, Pbar
  ENDTYPE T_TPCond
  
  TYPE:: T_TPSeries !container for a series of T,P Properties
    CHARACTER(15)::Name
    INTEGER:: NPoints
    REAL(dp),ALLOCATABLE:: vTdgC(:), vPbar(:)
  ENDTYPE T_TPSeries
  
  REAL(dp),PARAMETER:: Delta=1.0D0

CONTAINS

!!  IF(ALLOCATED(vTPpath)) DEALLOCATE(vTPpath)

INTEGER FUNCTION TPcond_IndexTP(TdgK,Pbar,vTPCond)
  USE M_Dtb_Const,ONLY: T_CK
  REAL(dp),       INTENT(IN):: TdgK,Pbar
  TYPE(T_TPCond), INTENT(IN):: vTPCond(:)
  INTEGER :: I,J
  J=0
  DO I=1,SIZE(vTPCond)
    IF(   ABS(TdgK -T_CK-vTPCond(I)%TdgC)<Delta &
    .AND. ABS(Pbar      -vTPCond(I)%Pbar)<Delta ) J=I
  ENDDO
  TPcond_IndexTP= J
ENDFUNCTION TPcond_IndexTP

INTEGER FUNCTION TPcond_IndexT(TdgK,vTPCond)
  USE M_Dtb_Const,ONLY: T_CK
  REAL(dp),      INTENT(IN):: TdgK
  TYPE(T_TPCond),INTENT(IN):: vTPCond(:)
  INTEGER :: I,J
  J=0
  DO I=1,SIZE(vTPCond)
    IF(ABS(TdgK -T_CK-vTPCond(I)%TdgC)<Delta) J=I
  ENDDO
  TPcond_IndexT= J
ENDFUNCTION TPcond_IndexT

ENDMODULE M_T_Tpcond

