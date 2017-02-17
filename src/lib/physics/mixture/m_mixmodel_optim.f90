module M_Mixmodel_Optim
!--
!-- interface between Mixture model and optimization routine
!--
  use M_Kinds
  use M_T_MixModel
  use M_Global_Vars
  !
  implicit none

  private

  !~ public:: Mixmodel_Optim_GetGibbs
  !~ public:: Mixmodel_Optim_GetPotentials
  public:: Mixmodel_Optim_GetMu
  public:: Mixmodel_Optim_SetParams
  !
  public:: Mixmodel_Optim_vMu0rt
  public:: Mixmodel_Optim_vLPole

  real(dp),allocatable:: Mixmodel_Optim_vMu0rt(:)
  logical, allocatable:: Mixmodel_Optim_vLPole(:)

  type(T_MixModel):: MixModel
  real(dp):: TdgK,Pbar

contains

subroutine Mixmodel_Optim_SetParams(T,P,MM) !,vG0rt)
  real(dp),intent(in):: T,P !T= dgK / P= bar
  type(T_MixModel),intent(in):: MM
  !real(dp),intent(in):: vG0rt(:)
  !
  TdgK= T
  Pbar= P
  MixModel= MM
  !
end subroutine Mixmodel_Optim_SetParams

!~ function Mixmodel_Optim_GetGibbs(vX) result(G)
  !~ real(dp),intent(in):: vX(:)
  !~ real(dp):: G
  !~ !
  !~ logical,allocatable:: vLPole(:)
  !~ !
  !~ real(dp):: Gmix,Gmeca
  !~ integer :: I
  !~ !
  !~ allocate(vLPole(size(vX)))
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
  !~ do I=1,MixModel%NPole
    !~ if(vLPole(I)) &
    !~ & Gmeca= Gmeca &
    !~ &      + vX(I) *vSpc(MixModel%vIPole(I))%G0rt
  !~ end do
  !~ !----------------------------------------------/"mechanical" mixing --
  !~ !
  !~ G= Gmix +Gmeca
  !~ !
  !~ deallocate(vLPole)
  !~ !
!~ end function Mixmodel_Optim_GetGibbs

subroutine Mixmodel_Optim_GetMu(vX,vMu)
  !real(dp),intent(in):: G
  real(dp),intent(in):: vX(:)
  real(dp),intent(out):: vMu(:)
  !
  integer:: N
  !
  real(dp),allocatable:: vLGam(:),vLIdeal(:),vLnAct(:)
  logical:: Ok
  character(len=30):: Msg
  !
  N= MixModel%NPole
  !
  allocate(vLGam(N),vLIdeal(N),vLnAct(N))
  !
  call MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MixModel,  & ! in: mixing model
  & vX,        & ! in
  & Mixmodel_Optim_vLPole, & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & !
  & vLIdeal,   & !
  & vLnAct)      !
  !
  where(Mixmodel_Optim_vLPole)  ;  vMu= vLnAct + Mixmodel_Optim_vMu0rt
  elsewhere                     ;  vMu= Zero
  end where
  !
  deallocate(vLGam,vLIdeal,vLnAct)
  !
  !do i=1,MixFas_Optim%Model%NPole
  !  vMu(I)= G0RT(I) +MixFas_Optim%vLnAct(I)
  !end do
  !
  !~ real(dp):: xx
  !~ real(dp),allocatable:: grad(:)
  !~ !
  !~ allocate(grad(size(vX)))
  !~ !
  !~ call gibbs_gradient(vX,grad)
  !~ xx= sum(vX(:)*grad(:))
  !~ vMu(:)= G - xx + grad(:)
  !~ !
  !~ deallocate(grad)
  !
end subroutine Mixmodel_Optim_GetMu

subroutine Mixmodel_Optim_GetPotentials(G,vX,vMu)
  real(dp),intent(in):: G
  real(dp),intent(in):: vX(:)
  real(dp),intent(out):: vMu(:)
  !
  integer:: i,N
  !
  logical, allocatable:: vLPole(:)
  real(dp),allocatable:: vLGam(:),vLIdeal(:),vLnAct(:)
  logical:: Ok
  character(len=30):: Msg
  !
  N= MixModel%NPole
  allocate(vLPole(N))
  vLPole(:)= (vX(:)>Zero) ! .and. (MixModel%vHasPole(:))
  !
  allocate(vLGam(N),vLIdeal(N),vLnAct(N))
  !
  call MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MixModel,  & ! in: mixing model
  & vX,        & ! in
  & vLPole,    & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & !
  & vLIdeal,   & !
  & vLnAct)      !
  !
  do i=1,N
    vMu(i)= vLnAct(i) + Mixmodel_Optim_vMu0rt(i)
  end do
  !
  deallocate(vLPole)
  deallocate(vLGam,vLIdeal,vLnAct)
  !
  !do i=1,MixFas_Optim%Model%NPole
  !  vMu(I)= G0RT(I) +MixFas_Optim%vLnAct(I)
  !end do
  !
  !~ real(dp):: xx
  !~ real(dp),allocatable:: grad(:)
  !~ !
  !~ allocate(grad(size(vX)))
  !~ !
  !~ call gibbs_gradient(vX,grad)
  !~ xx= sum(vX(:)*grad(:))
  !~ vMu(:)= G - xx + grad(:)
  !~ !
  !~ deallocate(grad)
  !
end subroutine Mixmodel_Optim_GetPotentials

end module M_Mixmodel_Optim
