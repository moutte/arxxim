module M_T_Kinmodel
!--
!-- data structure for kinetic model --
!--
  use M_Kinds
  
  implicit none
  
  private
  
  public:: MaxKinTerm
  public:: T_KinModel
  public:: KinModel_CalcCoeffs
  public:: KinModel_Index
  public:: KinModel_PrmIndex
  public:: KinModel_Show
  
  integer,parameter:: MaxKinTerm=4 !max nr of terms in the rate law
  
  type:: T_KinModel
    !implementation of a 'general' form of rate law, a la USGS (?)
    !rate/surf= SUM_i ( k_i * exp(-Ea_i /RT) * a_i^n_i * finh )  !with k_i=10^-pk_i
    !(finh=  inhibitor=  e.g. Al-species <-not implemented yet
    character(23):: Name !name of kinetic model
    !character(7):: Special !code for special case; not used now
    integer      :: NTermP,NTermD !number of terms, in Precip. and Dissol. laws, 
    !                             !up to MaxKinTerm each
    integer,      dimension(1:MaxKinTerm):: &
    & ISpcPrec, ISpcDiss, &    !index of "activator" species in current species array
    & JSpcPrec, JSpcDiss       !idem, for "double dependency", e.g. sulfides
    !
    real(dp),     dimension(1:MaxKinTerm)::&
    & Kd,pK_d,E_d,N_d,NJ_d,& !dissolution rate parameters !Kd_A=pKdA*function(-EdA/RT)
    & Kp,pK_p,E_p,N_p,NJ_p   !precipitation rate parameters !Kp_A=pKpA*function(-EdA/RT)
    !
    real(dp)::&
    & AlfaD,BetaD,& !VmQsK=(QsK**M%AlfaD - One)**M%BetaD, etc
    & AlfaP,BetaP
    !
  end type T_KinModel
  !
contains
  !
subroutine KinModel_CalcCoeffs(M,TimeFactor,TdgK) !-> Kd*exp(-Ea/RT)
!--
!-- temperature dependency of kinetic parameters --
!--
  use M_Dtb_Const,only: R_jK
  !
  type(T_KinModel),intent(inout)::M
  real(dp),        intent(in)   ::TimeFactor,TdgK
  !
  real(dp)::xRT
  !
  M%Kd=Zero
  M%Kp=Zero
  !
  ! PalandriKharaka, USGS_OpenFileReport, 2004-1068
  ! k(T)=k_298 * exp( - Ea_298 (1/T - 1/298) / R )
  ! M%Kd_A= TimeFactor * 10**(-M%pKdA) *exp(M%EdA*xRT)
  !
  xRT= (One/298.15D0 - One/TdgK) /R_jk
  !
  M%Kd(1:M%NTermD)= TimeFactor * 10**(-M%pK_d(1:M%NTermD)) *exp(M%E_d(1:M%NTermD)*xRT)
  M%Kp(1:M%NTermP)= TimeFactor * 10**(-M%pK_p(1:M%NTermP)) *exp(M%E_p(1:M%NTermP)*xRT)
  !
end subroutine KinModel_CalcCoeffs

integer function KinModel_Index(Str,V)
!--
!-- -> position of T_KinModel with %Name--Str in V
!--
  character(len=*),intent(in)::Str
  type(T_KinModel),intent(in)::V(:)
  !
  integer::I
  
  KinModel_Index=0
  !
  if(size(V)==0) return
  !
  I=0
  do
    I=I+1 
    if(trim(Str)==trim(V(I)%Name)) then
      KinModel_Index=I
      exit
    end if
    if(I==size(V)) exit
  end do ! if Str not found, KinModel_Index=0
  
  return
end function KinModel_Index

subroutine KinModel_PrmIndex(vPrmBk,M_In,M_Out)
!--
!-- update indexes of species in the kinetic model of phase M
!-- according to current species permutation
!--
  !
  integer,dimension(:),intent(in)   :: vPrmBk
  type(T_KinModel),    intent(in)   :: M_in
  type(T_KinModel),    intent(inout):: M_out
  !
  !-- precipitation law
  M_Out%ISpcPrec(1:M_Out%NTermP)= vPrmBk(M_In%ISpcPrec(1:M_In%NTermP))
  ! if another "activator" species
  where(M_Out%JSpcPrec(1:M_Out%NTermP) >0) &
  & M_Out%JSpcPrec(1:M_Out%NTermP)=vPrmBk(M_In%JSpcPrec(1:M_In%NTermP))
  !
  !-- dissolution law
  M_Out%ISpcDiss(1:M_Out%NTermD)=vPrmBk(M_In%ISpcDiss(1:M_In%NTermD))
  !-- when two activator species
  where(M_Out%JSpcDiss(1:M_Out%NTermD) >0) &
  & M_Out%JSpcDiss(1:M_Out%NTermD)=vPrmBk(M_In%JSpcDiss(1:M_In%NTermD))
  !
end subroutine KinModel_PrmIndex

subroutine KinModel_Show(f,vSpc,vKinMod,vPrm)
!--
  use M_T_Species,only: T_Species
  !
  integer,         intent(in):: F
  type(T_Species), intent(in):: vSpc(:)
  type(T_KinModel),intent(in):: vKinMod(:)
  integer,         intent(in),optional:: vPrm(:)
  !
  character:: T_=char(9)
  integer::I,J,K
  type(T_KinModel)::M
  !
  !call GetUnit(f)
  !open(f,file=trim(DirLog)//"kinmodel.log")
  !write(f,'(A,/)') "!.-> data on all minerals found in vSpc"
  !
  write(F,'(/,A)') "< KinModel_Show"
  write(F,'(A,/)')   "species, exponent, pK, ActivEnergy"
  
  do I=1,size(vKinMod)
    
    write(F,'(A)') trim(vKinMod(I)%Name)
    
    M=vKinMod(I)
    
    !-- dissolution terms --
    do J=1,M%NTermD
      
      K=M%ISpcDiss(J)
      
      if(present(vPrm)) K=vPrm(K) !for vSpc
      
      write(F,'(I3,2A,A1,3(F12.3,A1))') &
      & J,"_Dissol ",trim(vSpc(K)%NamSp),T_,M%N_D(J),T_,M%pK_D(J),T_,M%E_D(J),T_
      
      K=M%JSpcDiss(J)
      if(K>0) then
        if(present(vPrm)) K=vPrm(K) !for vSpc
        write(F,'(I3,2A,A1,F12.3)') &
        & J,"_Dissol ",trim(vSpc(K)%NamSp),T_,M%NJ_D(J) !,T_,M%pK_D(J),T_,M%E_D(J),T_
      end if
      
    end do
    
    !-- precipitation terms--
    do J=1,M%NTermP
    
      K=M%ISpcPrec(J)
      
      if(present(vPrm)) K=vPrm(K) !for vSpc
      
      write(F,'(I3,2A,A1,3(F12.3,A1))') &
      & J,"_Precip ",trim(vSpc(K)%NamSp),T_,M%N_P(J),T_,M%pK_P(J),T_,M%E_P(J),T_
      
      K=M%JSpcPrec(J)
      if(K>0) then
        if(present(vPrm)) K=vPrm(K) !for vSpc
        write(F,'(I3,2A,A1,3(F12.3,A1))') &
        & J,"_Precip ",trim(vSpc(K)%NamSp),T_,M%NJ_P(J) !,T_,M%pK_D(J),T_,M%E_D(J),T_
      end if
      
    end do
    
  end do
  write(F,'(A,/)') "</ KinModel_Show"
  
end subroutine KinModel_Show

end module M_T_Kinmodel 

