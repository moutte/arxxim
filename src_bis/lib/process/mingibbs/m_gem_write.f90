module M_GEM_Write
  !
  use M_Kinds
  use M_Trace, only: fTrc,iDebug,T_,Stop_
  !
  implicit none
  !
  private
  !
  public:: GEM_Write_Phases
  public:: GEM_Write_Mixtures
  public:: GEM_Write_Log_Entete
  !
contains

subroutine GEM_Write_Phases( &
& DimPath, &
& vFasIsPresent, &
& vSimplex_Ok, &
& tResult, &
& tResultMix)
!--
!-- write tabulated results for all phases found stable at any step
!--
  use M_Files,        only: DirOut
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: GetUnit
  use M_T_MixModel,   only: T_MixModel
  use M_GEM_Vars,     only: T_SavPhase
  !
  use M_Global_Vars,only: vSpc,vFas,vMixModel
  use M_GEM_Vars,   only: vCpnGEM
  !---------------------------------------------------------------------
  integer,       intent(in):: DimPath
  logical,       intent(in):: vFasIsPresent(:)
  logical,       intent(in):: vSimplex_Ok(DimPath)
  real(dp),      intent(in) :: tResult(:,:)
  type(T_SavPhase),intent(in):: tResultMix(:,:)
  !---------------------------------------------------------------------
  integer :: iPath,iFs,I,J,K
  integer :: nC,nF,nP,nFpur,nFmix
  integer :: F,F1,F2
  !real(dp):: vX(MaxPole)
  real(dp):: Tot,X,Y
  type(T_SavPhase):: Fas0
  type(T_MixModel):: MM
  !
  logical,parameter:: ComputeXMean= .true. !.false.
  !---------------------------------------------------------------------
  nC= size(vCpnGEM)
  nF= size(vFasIsPresent)
  nFpur= size(vFas)
  nFmix= size(vMixModel)
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_phase_mole.restab")
  call  Write_Title(F)
  !
  call GetUnit(F1)
  open(F1,file=trim(DirOut)//"_phase_cm3.restab")
  call  Write_Title(F1)
  !
  call GetUnit(F2)
  open(F2,file=trim(DirOut)//"_phase_grams.restab")
  call  Write_Title(F2)
  !
  K=0
  do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then
      !
      K=K+1
      !
      write(F, '(I3,A1)',      advance="no") K,T_  != count
      write(F1,'(I3,A1)',      advance="no") K,T_  != count
      write(F2,'(I3,A1)',      advance="no") K,T_  != count
      !
      write(F,'(2(F12.3,A1))',advance="no") &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      do iFs=1,nC != (3:nC+2)                     != components
        write(F,'(G15.6,A1)',advance="no") tResult(2+iFs,iPath),T_
      enddo
      !
      write(F1,'(2(F12.3,A1))',advance="no") &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      do iFs=1,nC != (3:nC+2)                     != components
        write(F1,'(G15.6,A1)',advance="no") tResult(2+iFs,iPath),T_
      enddo
      !
      write(F2,'(2(F12.3,A1))',advance="no") &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      do iFs=1,nC != (3:nC+2)                     != components
        write(F2,'(G15.6,A1)',advance="no") tResult(2+iFs,iPath),T_
      enddo
      !
      !~ Tot=SUM(tResult(1:nC,iPath))
      !~ do iFs=1,nC
        !~ write(F,'(G15.6,A1)',advance="no") tResult(iFs,iPath)/Tot,T_
      !~ enddo
      
      !-----------------------------------------------------mole numbers
      do iFs=1,nF != (nC+3:nC+nF+2)
        if(vFasIsPresent(iFs)) &
        & write(F,'(G15.3,A1)',advance="no") &
        ! nr'moles phase          *nr'oxygen in formula
        ! & tResult(nC+iFs,iPath)*tFormula(iFs,1),T_
        & tResult(2+nC+iFs,iPath),T_
      enddo
      write(F,*)
      !----------------------------------------------------/mole numbers

      !---------------------------------------------volume/cm3, weight/g
      do iFs=1,nFpur != (nC+3:nC+nF+2)
        if(vFasIsPresent(iFs)) &
        & write(F1,'(G15.3,A1)',advance="no") &
        & tResult(2+nC+iFs,iPath)*vFas(iFs)%VolM3*1.D6,T_
        if(vFasIsPresent(iFs)) &
        & write(F2,'(G15.3,A1)',advance="no") &
        & tResult(2+nC+iFs,iPath)*vFas(iFs)%WeitKg*1.D3,T_
      end do
      do iFs=1,nFmix
        if(vFasIsPresent(nFpur +iFs)) then
        
          Fas0= tResultMix(iFs,iPath)
          MM= vMixModel(Fas0%iModel)
          X= 0.D0
          Y= 0.D0
          do J=1,MM%nPole
            do I=1,MM%nPole
              X= X + vSpc(MM%vIPole(I))%V0     * Fas0%tXPole(J,I)
              Y= Y + vSpc(MM%vIPole(I))%WeitKg * Fas0%tXPole(J,I)
            end do
          end do
          
          !X= dot_product( &
          !& tResultMix(iFs,iPath)%tXPole(1,:), &
          !& tResultMix(iFs,iPath)%vVol0(:))
          
          write(F1,'(G15.3,A1)',advance="no") &
          & tResult(2+nC+nFpur+iFs,iPath)*X*1.D6,T_
          write(F2,'(G15.3,A1)',advance="no") &
          & tResult(2+nC+nFpur+iFs,iPath)*Y*1.D3,T_
          
        end if
      end do
      write(F1,*)
      write(F2,*)
      !--------------------------------------------/volume cm3, weight/g
    end if !!if(vSimplex_Ok(iPath))

  enddo !!do iPath
  !
  write(F,'(A1)') "_"
  !
  close(F)
  close(F1)
  close(F2)
  !
  if(iDebug>2) &
  & print '(2A)', "Phases: Results in ",trim(DirOut)//"_phase_xx.restab"
  !
contains

  subroutine Write_Title(FF)
    integer,intent(in):: FF
    !--------------------------------------------------------title lines
    write(FF,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
    !
    !--component names
    do i=1,nC
      write(FF,'(A,A1)',advance="no") trim(vCpnGEM(i)%NamCp)//"_cpn", T_
    enddo
    !
    !--pure phase names
    do iFs=1,nFpur
      if(vFasIsPresent(iFs)) &
      & write(FF,'(A,A1)',advance="no") trim(vFas(iFs)%NamFs), T_
    enddo
    !--mixture phase names
    do iFs=1,nFmix
      if(vFasIsPresent(nFpur +iFs)) &
      & write(FF,'(A,A1)',advance="no") trim(vMixModel(iFs)%Name), T_
    enddo
    write(FF,*)
    !-------------------------------------------------------/title lines
  end subroutine Write_Title

end subroutine GEM_Write_Phases

subroutine GEM_Write_Mixtures( &
& DimPath, &
& vFasIsPresent, &
& vSimplex_Ok,&
& TolX, &
& tResult, &
& tResultMix)
!--
!-- write tabulated results for all phases found stable at any step
!--
  use M_Files,        only: DirOut
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: GetUnit
  use M_Numeric_Tools,only: iMaxLoc_R
  !
  use M_T_Phase,      only: T_Phase
  use M_T_MixModel,   only: T_MixModel,MaxPole,Mix_Molecular
  !
  use M_Global_Vars,  only: vSpc,vFas,vMixModel
  use M_GEM_Vars,     only: T_SavPhase
  use M_GEM_Vars,     only: vCpnGEM
  !---------------------------------------------------------------------
  integer, intent(in):: DimPath
  logical, intent(in):: vFasIsPresent(:)
  logical, intent(in):: vSimplex_Ok(DimPath)
  real(dp),intent(in):: TolX ! <-MixMinim_TolX
  real(dp),intent(in) :: tResult(:,:) ! for T,P values
  type(T_SavPhase),intent(inout):: tResultMix(:,:)
  !---------------------------------------------------------------------
  integer :: iPath,iFs,I,J,K,P,Q
  integer :: nP,nFmix,nFPur,nPP
  integer :: FMIX
  real(dp):: vX(MaxPole)
  real(dp):: Tot
  type(T_SavPhase):: Fas0,Fas1,SavPhaseZero
  type(T_MixModel):: MM
  !
  logical,parameter:: ComputeXMean= .true. !.false.
  !---------------------------------------------------------------------
  call GetUnit(FMIX)
  open(FMIX,file=trim(DirOut)//"_mixtures.restab")
  !
  nFpur= size(vFas)
  nFmix= size(vMixModel)
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  !-----------------------------------------------------------title line
  write(FMIX,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
  !
  do iFs=1,nFmix
    !
    if(vFasIsPresent(nFpur +iFs)) then
      !
      !Fas= tResultMix(iFs,iPath)
      MM= vMixModel(iFs)
      nP= MM%nPole
      !if(Fas0%iModel /= 0) then
      if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
        nPP= 1
      else
        nPP= MM%nPole
      end if
      do I=1,nPP
        write(FMIX,'(A,A1)',advance="NO") trim(MM%Name),T_
        do P=1,nP
          write(FMIX,'(A,A1)',advance="NO") &
          & trim(vSpc(MM%vIPole(P))%NamSp),T_
        enddo
      enddo
      !
    end if
  
  enddo
  write(FMIX,*)
  !----------------------------------------------------------/title line
  
  K=0
  do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then
      !
      K=K+1
      !------------------------------compute mean comp'n of each mixture
      !----------------------------------------------and its free energy
      if(ComputeXMean) then
        !
        do J=1,nFmix
          !
          Fas0=        tResultMix(J,iPath)
          Fas1=        SavPhaseZero
          Fas1%iModel= Fas0%iModel
          nP=          vMixModel(Fas0%iModel)%nPole
          Fas1%nFas=   nP
          !
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          !if(Fas0%iModel /= 0) then
          do P=1,Fas0%nFas
            do Q=1,Fas1%nFas

              vX(1:nP)= ABS( Fas0%tXPole(P,1:nP) - Fas0%tXPole(Q,1:nP) )

              if(MAXVAL(vX(1:nP)) < TolX *1.0D2) then
              
                Tot= Fas0%vMole(P) +Fas1%vMole(Q)
                vX(1:nP)= Fas0%vMole(P) *Fas0%tXPole(P,1:nP) &
                &       + Fas1%vMole(Q) *Fas0%tXPole(Q,1:nP)
                Fas1%tXPole(Q,1:nP)= vX(1:nP) /Tot
                Fas1%vMole(Q)= Tot
                
                exit
              
              end if

            enddo
          enddo
          !
          tResultMix(J,iPath)= Fas1
          !end if
          !
        enddo
        !
      end if
      !----------------------------------/compute mean comp'n of mixture

      !---------------------------------------write mixture compositions
      !-----------------------------------------------------in file FMIX
      write(FMIX,'(I3,A1)',     advance="no") K,T_   != count
      !
      write(FMIX,'(2(F7.2,A1))',advance="no") &      !
      & tResult(1,iPath), T_, &                      != TdgC
      & tResult(2,iPath), T_                         != Pbar

      do J=1,nFmix
        if(vFasIsPresent(nFpur +J)) then
          Fas0= tResultMix(J,iPath)
          MM= vMixModel(Fas0%iModel)
          nP= MM%nPole
          !if(Fas0%iModel /= 0) then
          if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
            nPP= 1
          else
            nPP= MM%nPole
          end if
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          do I=1,nPP
            write(FMIX,'(G15.6,A1)',advance="NO") Fas0%vMole(I),T_
            do P=1,nP
              write(FMIX,'(F7.3,A1)',advance="NO") &
              & Fas0%tXPole(I,P),T_
            end do
          end do
          !end if
        end if
      enddo
      !
      write(FMIX,*)
      !--------------------------------------/write mixture compositions
      !
    end if !!if(vSimplex_Ok(iPath))

  enddo !!do iPath
  !
  if(iDebug>2) &
  & print '(2A)', "Mixtures: Results in ",trim(DirOut)//"_mixtures.restab"
  !
end subroutine GEM_Write_Mixtures

subroutine GEM_Write_Log_Entete(vFas)
  use M_IoTools,only: GetUnit
  use M_Files,  only: DirOut
  use M_T_Phase,only: T_Phase
  !
  type(T_Phase),intent(in):: vFas(:)
  !
  integer:: F1,F2
  !
  call GetUnit(F1)
  open(F1,file=trim(DirOut)//"_gibbs.log")
  write(F1,'(A)') "g of all phases, relative to current assemblage"
  !
  call GetUnit(F2)
  open(F2,file=trim(DirOut)//"_mixcomp.log")
  write(F2,'(A)')   "properties of mixtures with min(g)<0 at each iteration"
  write(F2,'(A,/)') "mix.model.name, Xi, Gmin"
  !
end subroutine GEM_Write_Log_Entete

subroutine Simplex_WriteTable(iFil,M,N)
  use M_IoTools,     only: OutStrVec 
  use M_Simplex_Vars,only: tSimplex 
  integer,intent(in)::iFil,M,N
  !
  real(dp),dimension(1:N+1)::V
  integer:: iLin !, iCol
  !
  do iLin=0,M+1
    V(1:N+1)=tSimplex(iLin,0:N)
    call OutStrVec(iFil,V,Opt_C="7")
  enddo
  !
  write(iFil,'(A1)') "_"
  !
end subroutine Simplex_WriteTable

end module M_GEM_Write

!Simplex, notes from NumRecp:
!
!linear programming= 
!maximize the linear function 
!  Z=  a(0,1).x(1) + ... + a(0,n).x(n)
!subject to 
!  primary constraints: a(0,1:n)>=0
!and simultaneously subject to 
!  M=  m1 + m2 + m3 additional constraints,
!  m1 of them of the form: a(i,1).x1 + ... + a(i,N).x(N) >= b(i)   //i=1..m1
!  m2 of them of the form: a(j,1).x1 + ... + a(j,N).x(N) <= b(j)   //j=m1+1..m1+m2
!  m3 of them of the form: a(k,1).x1 + ... + a(k,N).x(N) <= b(k)   //k=m1+m2+1..m1+m2+m3

!On output, the tableau a is indexed by two returned arrays of integers
!iposv(j) contains, for j= 1...M, 
!the number i whose original variable x(i) is now represented by row j+1 of a
!These are thus the left-hand variables in the solution
!(The first row of a is of course the z-row.)
!A value i > N indicates that the variable is a y(i) rather than an x(i), x(N+j)==y(j)
!Likewise, izrov(j) contains, for j= 1..N,
!the number i whose original variable x(i) is now a right-hand variable, represented by column j+1 of a. 
!These variables are all zero in the solution. 
!The meaning ofi > N is the same as above, 
!except that i >N +m1 +m2 denotes an artificial or slack variable 
!which was used only internally and should now be entirely ignored.


