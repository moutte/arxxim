module M_Dynam_Solve_Vars
!--
!-- dynamic variables used only at "low level" (residual / jacobian)
!-- -> separated from other variable sets for dynam modules
!-- (M_Dynam_Vars,M_Dynam_OneStep_Vars)
!--
  use M_Kinds
  implicit none

  private
  !
  public:: Dynam_Solve_Vars_Alloc
  public:: Dynam_Solve_Vars_Clean
  !
  real(dp),dimension(:,:),allocatable,public:: &
  !! & dVmAdLnXf, &   !1:nMk,1:nAq, @VmAct / @LnXf
  & dVmQdLnXf, &   !1:nMk,1:nCi, @VmQsK / @LnXf
  ! 
  & dSRdLnXm,  &   !1:nMk,1:nMk, derivative of vSurfK / Ln(nMol(Mineral))
  & dSRdXm,    &   !1:nMk,1:nMk, derivative of vSurfK / nMol(Mineral)
  !
  & dVMdLnXf,  &   !1:nMk,1:nCi, derivative of vRateM / Ln(nMol(Prim'Aqu.Spcecies))
  & dVMdLnXm,  &   !1:nMk,1:nMk, derivative of vRateM / Ln(nMol(Mineral))
  & dVMdXm         !1:nMk,1:nMk, derivative of vRateM / nMol(Mineral)
  !
contains

subroutine Dynam_Solve_Vars_Alloc(nMk,nCi)
  integer, intent(in):: nMk,nCi
  !
  allocate(dSRdLnXm(1:nMk,1:nMk))     ;  dSRdLnXm= Zero
  allocate(dSRdXm(1:nMk,1:nMk))       ;  dSRdXm=   Zero
  !
  allocate(dVmQdLnXf(1:nMk,1:nCi))    ;  dVmQdLnXf= Zero !VmAct depends on fluid species only
  !
  !! allocate(dVmAdLnXf(1:nMk,1:nCi)); dVmAdLnXf= Zero !VmAct depends on fluid species only
  !
  allocate(dVmdLnXf(1:nMk,1:nCi))     ;  dVmdLnXf= Zero !
  allocate(dVmdLnXm(1:nMk,1:nMk))     ;  dVmdLnXm= Zero !
  allocate(dVmdXm  (1:nMk,1:nMk))     ;  dVmdXm=   Zero !
  !
end subroutine Dynam_Solve_Vars_Alloc

subroutine Dynam_Solve_Vars_Clean

  if(allocated(dSRdLnXm))   deallocate(dSRdLnXm)
  if(allocated(dSRdXm))     deallocate(dSRdXm)
  !
  !! if(allocated(dVmAdLnXf))  deallocate(dVmAdLnXf)      
  if(allocated(dVmQdLnXf))  deallocate(dVmQdLnXf)
  !
  if(allocated(dVmdXm))     deallocate(dVmdXm)
  if(allocated(dVmdLnXf))   deallocate(dVmdLnXf) 
  if(allocated(dVmdLnXm))   deallocate(dVmdLnXm)

end subroutine Dynam_Solve_Vars_Clean

end module M_Dynam_Solve_Vars
