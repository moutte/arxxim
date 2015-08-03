MODULE M_Mixmodel_Optim
!--
!-- interface between Mixture model and optimization routine
!--
  USE M_Kinds
  USE M_T_MixModel
  USE M_Global_Vars
  !
  IMPLICIT NONE

  PRIVATE

  !~ PUBLIC:: Mixmodel_Optim_GetGibbs
  !~ PUBLIC:: Mixmodel_Optim_GetPotentials
  PUBLIC:: Mixmodel_Optim_GetMu
  PUBLIC:: Mixmodel_Optim_SetParams
  !
  PUBLIC:: Mixmodel_Optim_vMu0rt
  PUBLIC:: Mixmodel_Optim_vLPole

  REAL(dp),ALLOCATABLE:: Mixmodel_Optim_vMu0rt(:)
  LOGICAL, ALLOCATABLE:: Mixmodel_Optim_vLPole(:)

  TYPE(T_MixModel):: MixModel
  REAL(dp):: TdgK,Pbar

CONTAINS

SUBROUTINE Mixmodel_Optim_SetParams(T,P,MM) !,vG0rt)
  REAL(dp),INTENT(IN):: T,P !T= dgK / P= bar
  TYPE(T_MixModel),INTENT(IN):: MM
  !REAL(dp),INTENT(IN):: vG0rt(:)
  !
  TdgK= T
  Pbar= P
  MixModel= MM
  !
END SUBROUTINE Mixmodel_Optim_SetParams

!~ FUNCTION Mixmodel_Optim_GetGibbs(vX) RESULT(G)
  !~ REAL(dp),INTENT(IN):: vX(:)
  !~ REAL(dp):: G
  !~ !
  !~ LOGICAL,ALLOCATABLE:: vLPole(:)
  !~ !
  !~ REAL(dp):: Gmix,Gmeca
  !~ INTEGER :: I
  !~ !
  !~ ALLOCATE(vLPole(SIZE(vX)))
  !~ vLPole(:)= (vX(:)>Zero)
  !~ !
  !~ Gmix= MixModel_GibbsMixRT( & !
  !~ & TdgK,Pbar, & !IN
  !~ & MixModel,  & !IN, mixing model
  !~ & vLPole,    & !IN
  !~ & vX)          !IN
  !~ !
  !~ !---------------------------------------------- "mechanical" mixing --
  !~ Gmeca= Zero
  !~ DO I=1,MixModel%NPole
    !~ IF(vLPole(I)) &
    !~ & Gmeca= Gmeca &
    !~ &      + vX(I) *vSpc(MixModel%vIPole(I))%G0rt
  !~ ENDDO
  !~ !----------------------------------------------/"mechanical" mixing --
  !~ !
  !~ G= Gmix +Gmeca
  !~ !
  !~ DEALLOCATE(vLPole)
  !~ !
!~ END FUNCTION Mixmodel_Optim_GetGibbs

SUBROUTINE Mixmodel_Optim_GetMu(vX,vMu)
  !REAL(dp),INTENT(IN):: G
  REAL(dp),INTENT(IN):: vX(:)
  REAL(dp),INTENT(OUT):: vMu(:)
  !
  INTEGER:: N
  !
  REAL(dp),ALLOCATABLE:: vLGam(:),vLIdeal(:),vLnAct(:)
  LOGICAL:: Ok
  CHARACTER(LEN=30):: Msg
  !
  N= MixModel%NPole
  !
  ALLOCATE(vLGam(N),vLIdeal(N),vLnAct(N))
  !
  CALL MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MixModel,  & ! in: mixing model
  & vX,        & ! in
  & Mixmodel_Optim_vLPole, & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & !
  & vLIdeal,   & !
  & vLnAct)      !
  !
  WHERE(Mixmodel_Optim_vLPole)  ;  vMu= vLnAct + Mixmodel_Optim_vMu0rt
  ELSEWHERE                     ;  vMu= Zero
  END WHERE
  !
  DEALLOCATE(vLGam,vLIdeal,vLnAct)
  !
  !DO i=1,MixFas_Optim%Model%NPole
  !  vMu(I)= G0RT(I) +MixFas_Optim%vLnAct(I)
  !END DO
  !
  !~ REAL(dp):: xx
  !~ REAL(dp),ALLOCATABLE:: grad(:)
  !~ !
  !~ allocate(grad(size(vX)))
  !~ !
  !~ call gibbs_gradient(vX,grad)
  !~ xx= sum(vX(:)*grad(:))
  !~ vMu(:)= G - xx + grad(:)
  !~ !
  !~ deallocate(grad)
  !
END SUBROUTINE Mixmodel_Optim_GetMu

SUBROUTINE Mixmodel_Optim_GetPotentials(G,vX,vMu)
  REAL(dp),INTENT(IN):: G
  REAL(dp),INTENT(IN):: vX(:)
  REAL(dp),INTENT(OUT):: vMu(:)
  !
  INTEGER:: i,N
  !
  LOGICAL, ALLOCATABLE:: vLPole(:)
  REAL(dp),ALLOCATABLE:: vLGam(:),vLIdeal(:),vLnAct(:)
  LOGICAL:: Ok
  CHARACTER(LEN=30):: Msg
  !
  N= MixModel%NPole
  ALLOCATE(vLPole(N))
  vLPole(:)= (vX(:)>Zero) ! .AND. (MixModel%vHasPole(:))
  !
  ALLOCATE(vLGam(N),vLIdeal(N),vLnAct(N))
  !
  CALL MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MixModel,  & ! in: mixing model
  & vX,        & ! in
  & vLPole,    & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & !
  & vLIdeal,   & !
  & vLnAct)      !
  !
  DO i=1,N
    vMu(i)= vLnAct(i) + Mixmodel_Optim_vMu0rt(i)
  ENDDO
  !
  DEALLOCATE(vLPole)
  DEALLOCATE(vLGam,vLIdeal,vLnAct)
  !
  !DO i=1,MixFas_Optim%Model%NPole
  !  vMu(I)= G0RT(I) +MixFas_Optim%vLnAct(I)
  !END DO
  !
  !~ REAL(dp):: xx
  !~ REAL(dp),ALLOCATABLE:: grad(:)
  !~ !
  !~ allocate(grad(size(vX)))
  !~ !
  !~ call gibbs_gradient(vX,grad)
  !~ xx= sum(vX(:)*grad(:))
  !~ vMu(:)= G - xx + grad(:)
  !~ !
  !~ deallocate(grad)
  !
END SUBROUTINE Mixmodel_Optim_GetPotentials

END MODULE M_Mixmodel_Optim
