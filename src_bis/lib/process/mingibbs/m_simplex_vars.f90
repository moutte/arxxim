module M_Simplex_Vars
!--
!-- vars for simplex computations
!--
  use M_Kinds
  use M_T_MixModel,only: MaxPole
  implicit none
  !
  private
  !
  public:: Simplex_Vars_Alloc
  public:: Simplex_Vars_Clean
  public:: tSimplex,IZROV,IPOSV
  !
  real(dp),allocatable:: tSimplex(:,:)
  integer, allocatable:: IZROV(:),IPOSV(:)

contains

subroutine Simplex_Vars_Alloc(M,N)
  integer,intent(in)::M,N
  !
  call Simplex_Vars_Clean
  !
  allocate(tSimplex(0:M+1,0:N))  ;  tSimplex=Zero
  allocate(IPOSV(1:M))           ;  IPOSV= 0
  allocate(IZROV(1:N))           ;  IZROV= 0
  !
end subroutine Simplex_Vars_Alloc

subroutine Simplex_Vars_Clean
  if(allocated(tSimplex))  deallocate(tSimplex)
  if(allocated(IZROV))     deallocate(IZROV)
  if(allocated(IPOSV))     deallocate(IPOSV)
end subroutine Simplex_Vars_Clean

end module M_Simplex_Vars
