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

subroutine SolModel_Alloc
  use M_T_SolModel
  use M_Global_Vars,only: vSolModel
  !
  integer:: nSolModel
  integer:: I
  
  nSolModel= 1
  
  if(allocated(vSolModel)) deallocate(vSolModel)
  allocate(vSolModel(nSolModel))
  
  do I=1,nSolModel
    !~ vSolModel(I)%Name= trim(vSolModelName(I))
    vSolModel(I)%Typ=  "LIQ"
    !~ vSolModel(I)%ActModel= "IDEAL"
    vSolModel(I)%iActModel= 1
  enddo
  
end subroutine SolModel_Alloc

subroutine SolPhase_Alloc
  use M_T_SolPhase
  use M_Global_Vars,only: vSolFas
  !
  if(allocated(vSolFas)) deallocate(vSolFas)
  allocate(vSolFas(1))
  
  vSolFas(1)%Name=   "WATER"
  vSolFas(1)%iModel= 1
  
end subroutine SolPhase_Alloc

end module M_SolModel_Alloc
