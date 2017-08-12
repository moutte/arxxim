module M_SolModel_Alloc
  use M_Kinds
  !
  implicit none
  !
  private
  !
  public:: SolModel_Alloc
  public:: SolPhase_Alloc
  !
contains

subroutine SolModel_Alloc(n,vSolModel)
  use M_T_SolModel
  ! use M_Global_Vars,only: vSolModel
  !
  integer,                     intent(in):: n
  type(T_SolModel),allocatable,intent(out):: vSolModel(:)
  !
  integer:: I
  
  ! if(allocated(vSolModel)) deallocate(vSolModel)
  allocate(vSolModel(n))
  
  do I=1,n
    !! vSolModel(I)%Name= trim(vSolModelName(I))
    vSolModel(I)%Typ=  "LIQ"
    !! vSolModel(I)%ActModel= "IDEAL"
    vSolModel(I)%iActModel= 1
  end do
  
end subroutine SolModel_Alloc

subroutine SolPhase_Alloc(n,vSolFas)
  use M_T_SolPhase
  !
  integer,                     intent(in):: n
  type(T_SolPhase),allocatable,intent(out):: vSolFas(:)
  !
  allocate(vSolFas(n))
  
end subroutine SolPhase_Alloc

end module M_SolModel_Alloc
