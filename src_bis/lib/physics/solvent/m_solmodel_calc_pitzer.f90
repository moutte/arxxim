module M_Solmodel_Calc_Pitzer
!--
!-- Module routines associées à Pitzer
!-- Livraison 10/2006: Nicolas Ferrando
!-- Remoduling : Anthony Michel
!-- Modified : J.Moutte, 11/2006, 01/2014
!-- -> M_Solmodel_Pitzer_Dtb-  database reading, build vPitz, vIPitz
!-- -> M_Solmodel_Pitzer_Calc- calculations
!--

  use M_Kinds
  use M_Trace,only: fTrc,iDebug,T_
  
  implicit none
  
  private
  !
  public:: Solmodel_Calc_Pitzer
  public:: EThetaCalc
  
  !! public:: Pitzer_Calc_APhiMonnin
  
  integer:: fTrcPitz
  
contains

subroutine Solmodel_Calc_Pitzer( &
& vSpc,                     & !IN
& MWsv,Rho,Eps,TdgK,vMolal, & !IN
& vLnGam,LnActSv,Osmo)        !OUT

  use M_Numeric_Const,only: Pi
  use M_T_Species,    only: T_Species
  !
  type(T_Species),intent(in) :: vSpc(:)
  real(dp),intent(in) :: MWsv,Rho,Eps,TdgK
  real(dp),intent(in) :: vMolal(:)
  !
  real(dp),intent(out):: vLnGam(:)
  real(dp),intent(out):: LnActSv   !solvent activity
  real(dp),intent(out):: Osmo
  !
  type(T_Species),allocatable:: vSolut(:)
  real(dp),       allocatable:: vMolalSolut(:)
  integer :: I,J,N
  real(dp):: APhi
  
  call Pitzer_Calc_APhi( &
  & Rho,Eps,TdgK, &
  & APhi)
  ! result at 25C/1ATM = 0.392521
  !
  ! call Pitzer_Calc_APhiMonnin(TdgK,1.013D0,Aphi)
  ! result at 25C/1ATM = 0.3914566
  !  
  ! APhi= 2.303D0 *0.5092D0 / 3.D0 !! Aphi
  ! = 0.380896
  !
  !! comparison of values of APhi at 298K/1atm
  !! 1- according to Pitzer formula    = 0.392521
  !! 2- according to Monnin regression = 0.391457
  !! 3- according to Bethke            = 0.380896
  
  N= count(vSpc(:)%Typ=="AQU") !- 1 !mod 19/06/2008 09:38
  allocate(vSolut(N))
  allocate(vMolalSolut(N))
  
  J= 0
  do I=1,size(vSpc)
    if(vSpc(I)%Typ=="AQU" .and. vSpc(I)%NamSp/="H2O") then
      J= J+1
      vSolut(J)= vSpc(I)
      vMolalSolut(J)= vMolal(I)
    end if
  enddo
  
  call Pitzer_Calc_Gamma( &
  & vSolut,Aphi,vMolalSolut, & !IN
  & Osmo,vLnGam) !OUT
  !
  LnActSv= -Osmo *MWsv * (SUM(vMolalSolut))
  
  deallocate(vSolut)
  deallocate(vMolalSolut)
  
  return
end subroutine Solmodel_Calc_Pitzer

subroutine Pitzer_Calc_Gamma( &
& vSpc,Aphi,vMolal, &
& Osmo,vLnGam) !,ERR)
!--
!-- Calcule les coefficients d'activité des solute's
!-- et le coefficient osmotique de l'eau par le mod'ele de Pitzer
!--
!-- vSpc is the list of solute species only
!--
!--IN
!--  vSpc:   solute's
!--  Aphi:   pente du terme de Debye-Huckel
!--  vMolal: molalite's des solute's
!--
!--OUT
!--  Osmo: Coefficient osmotique 
!--  vLnGam: log neperien des coefficient d'activite des constituants
!--  ERR: indicateur d'erreur (non utilise -> a completer si besoin)
!--
  use M_IoTools,only:GetUnit
  use M_T_Species,only: T_Species
  use M_Solmodel_Pitzer_Dtb,only: tBeta0,tBeta1,tBeta2,tAlfa1,tAlfa2
  use M_Solmodel_Pitzer_Dtb,only: tTheta,tPsi,tCPhi,tZeta,tLamda
  
  type(T_Species),intent(in) :: vSpc(:)
  real(dp),       intent(in) :: Aphi
  real(dp),       intent(in) :: vMolal(:)
  !
  real(dp),       intent(out):: Osmo
  real(dp),       intent(out):: vLnGam(:)
  !
  !------------------------------------------------ Variables locales --
  integer :: iC,jC,kC,iA,jA,kA,jN,A1,A2,C1,C2,K1,K2
  integer :: I,J,K,nSp
  integer :: FF !file index
  real(dp):: Ionic, M, BigZ, BigF
  real(dp):: RI, Fi0, Fi1, fOsmo
  real(dp):: AC, Tmp, X
  real(dp),dimension(1:size(vSpc),1:size(vSpc)):: &
  & Bmx1, & ! B_MX  of Pitzer equations, p117, Bethke'96
  & Bmx2, & ! B'_MX of Pitzer equations, p117, Bethke'96
  & Cmx,  & ! C_MX  of Pitzer equations, p118, Bethke'96
  & tETheta, tETheta1
  !
  integer,allocatable:: vCharge(:),vIonType(:)
  !---------------------------------------------------------------------
  
  if(iDebug>0) write(fTrc,'(A)') "< Pitzer_Calc_Gamma"
  !
  nSp= size(vSpc)
  !
  allocate(vCharge(1:nSp))
  vCharge(:)=vSpc(:)%Z
  !
  allocate(vIonType(1:nSp))
  vIonType(:)= 0
  where(vCharge(:)>0) vIonType(:)=  1
  where(vCharge(:)<0) vIonType(:)= -1
  !
  if(iDebug>2) then
    !
    call GetUnit(ff)
    open(ff,file="debug_pitzer.log")
    !
    write(ff,'(A,/)') "trace of routine Pitzer_Calc_Gamma"
    !
    call CheckCoeffs
    !
  end if
  
  !---------------------------------------------------- FORCE IONIQUE --
  M=    Zero  !! total molality
  BigZ= Zero  !! total molal charge, sum(|z_i|.m_i)
  Ionic=Zero  !! ionic strength
  do I=1,nSp
    !print *,I
    M=     M     + vMolal(I)
    BigZ=  BigZ  + vMolal(I)*vCharge(I)*vIonType(I)
    Ionic= Ionic + vMolal(I)*vCharge(I)*vCharge(I)
    !print *,Ionic
  enddo
  !pause
  Ionic= Ionic/2.0D0
  RI=  SQRT(Ionic)
  !
  if(iDebug>2) write(ff,'(A,G15.6)') "Ionic= ", Ionic
  !---------------------------------------------------/ FORCE IONIQUE --
  
  !--------------------------------------- COEFFICIENTS LONGUE PORTEE --
  Fi0= -4.0d0*Aphi*Ionic/1.2D0*log(One + 1.2D0*RI)
  Fi1= -4.0d0*Aphi      /1.2D0*log(One + 1.2D0*RI) &
  &    -2.0d0*Aphi*RI/(One + 1.2D0*RI)
  if(iDebug>2) write(ff,'(A,G15.6)') "Fi0=   ", Fi0
  if(iDebug>2) write(ff,'(A,G15.6)') "Fi1=   ", Fi1
  !--------------------------------------/ COEFFICIENTS LONGUE PORTEE --
  
  !------------------------------------------------ tETheta, tETheta1 --
  tETheta(:,:)= Zero
  tETheta1(:,:)= Zero
  do I=1,nSp
    do J=1,nSp
      !
      if(I==J) cycle
      if(vIonType(I)*vIonType(J)==1) then ! ions same sign
        !
        call EThetaCalc( &
        & vCharge(I),vCharge(J), RI, Aphi, &
        & tETheta(I,J), tETheta1(I,J))
        !
        if(iDebug>2) write(ff,'(A,2(G15.6,1X),2(A,A1))') &
        & "ETheta, ETheta1",&
        & tETheta(I,J), tETheta1(I,J), &
        & trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_
        !
      end if
      !
    end do
  end do
  !------------------------------------------------/tETheta, tETheta1 --
  
  Bmx1= Zero
  Bmx2= Zero
  Cmx=  Zero
  !
  !-------------------------------------- PARAMETRES DE PITZER B et C --
  !------------------------------------------------------ B_MX, B'_MX --
  !------------------------------ using tBeta0, tBeta1, tBeta2, tCPhi --
  if (Ionic /= Zero) then
    do I=1,nSp
      do J=1,nSp
        if(I==J) cycle
        !
        !X1=tAlfa1(I,J)*RI
        !X2=tAlfa2(I,J)*RI
        !
        Bmx1(I,J)= tBeta0(I,J) &
        &        + tBeta1(I,J)*FONC1(tAlfa1(I,J)*RI) & 
        &        + tBeta2(I,J)*FONC1(tAlfa2(I,J)*RI)
        !
        Bmx2(I,J)=(tBeta1(I,J)*FONC2(tAlfa1(I,J)*RI)    &
        &        + tBeta2(I,J)*FONC2(tAlfa2(I,J)*RI)  ) &
        &        /Ionic
        !
        if (vCharge(I)*vCharge(J)/=0) & !C_cation-anion eq24, Moeller'98
        & Cmx(I,J)=  0.5D0*tCPhi(I,J) / SQRT(ABS(real(vCharge(I)*vCharge(J))))
        !
      end do
    end do
    !
    !------------------------------------------------------------trace--
    if(iDebug>2) then
      do I=1,nSp
        do J=1,nSp
          if(I==J) cycle
          if(Bmx1(I,J)/=Zero) &
          & write(ff,'(A,3(G15.6,1X),2(A,A1))') &
          & "Beta-0-1-2   ", &
          & tBeta0(I,J), tBeta1(I,J), tBeta2(I,J), &
          & trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_
          !
        end do
      end do
      do I=1,nSp
        do J=1,nSp
          if(I==J) cycle
          if(Bmx1(I,J)/=Zero) &
          & write(ff,'(A,3(G15.6,1X),2(A,A1))') &
          & "Bmx1,Bmx2,Cmx", &
          & Bmx1(I,J), Bmx2(I,J), Cmx(I,J), &
          & trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_
        end do
      end do
    end if
  end if
  !---/
  !--------------------------------------/PARAMETRES DE PITZER B et C --
  
  !---------------------------------------------------------------------
  !-------------------------------------------- COEFFICIENT OSMOTIQUE --
  if (Ionic>Zero) then
  
    fOsmo= -Aphi*RI**3 /(One + 1.2D0*RI)
    
    if(iDebug>2) write(ff,'(I3,A,G15.6)') I," fOsmo=", fOsmo
    !fOsmo=  Ionic*Fi1 - Fi0
    !
    !------------------------------------- somme croisee anion/cation --
    tmp= Zero
    do I=1,nSp
      do J=1,nSp
        if(I==J) cycle
        if(vIonType(I)*vIonType(J)==-1) then ! anion/cation
          tmp= tmp &
          &  + vMolal(I)*vMolal(J) &
          &    * ( Bmx1(I,J) + Ionic*Bmx2(I,J) + BigZ*Cmx(I,J) )
        end if
      enddo
    enddo
    if(iDebug>2) write(ff,'(I3,A,G15.6)') I," fOsmo, anion/cation=", tmp
    fOsmo= fOsmo + tmp
    !call Pause_
    
    !---------------------------------- somme croisee Cation1/Cation2 --
    !--------------------------------------------- using tPsi, tTheta --
    tmp= Zero
    do C1=1,nSp-1
      if(vIonType(C1)/=1) cycle     ! C1 is CATION
      do C2=C1+1,nSp
        if(vIonType(C2)/=1) cycle   ! C2 is CATION
        !
        ! SUM_ON_ANIONS PSI(C1,C2,ANI) * vMolal(ANI)
        X=  Zero
        do K=1,nSp
          if(vIonType(K)==-1) &     !K is ANION
          & X= X + vMolal(K) *tPsi(C1,C2,K)
        enddo
        !X=  Zero
        !
        tmp= tmp &
        &  + vMolal(C1)*vMolal(C2) &
        &    *(tTheta(C1,C2) +tETheta(C1,C2) +RI*RI*tETheta1(C1,C2) +X)
        !
      enddo
      !
    enddo
    if(iDebug>2) write(ff,'(I3,A,G15.6)') C1," fOsmo, cation1/cation2=", tmp
    fOsmo= fOsmo + tmp
    !-----------------------------------------------/ Cation1/Cation2 --
    
    !------------------------------------ Somme Croisee Anion1/Anion2 --
    !--------------------------------------------- using tPsi, tTheta --
    tmp= Zero
    do A1=1,nSp-1
      if(vIonType(A1)/=-1) cycle     ! A1 is ANION
      do A2=A1+1,nSp
        if(vIonType(A2)/=-1) cycle   ! A2 is ANION
        !
        ! SUM_ON_CATIONS PSI(A1,A2,CAT) * vMolal(CAT)
        X=  Zero
        do K=1,nSp
          if(vIonType(K)==1) & !K is cation
          & X= X + vMolal(K) *tPsi(A1,A2,K)
        enddo
        ! X=  Zero
        !
        ! BigPhi^phi_ij= theta_ij + ETheta_ij + IStrength * EThetaPrim_ij
        tmp= tmp &
        &    + vMolal(A1) *vMolal(A2) &
        &      *(tTheta(A1,A2) +tETheta(A1,A2) +RI*RI*tETheta1(A1,A2) +X)
        !
      enddo
      !
    enddo
    if(iDebug>2) write(ff,'(I3,A,G15.6)') A1," fOsmo, Anion1/Anion2=", tmp
    fOsmo= fOsmo + tmp
    !--------------------------------------------------/Anion1/Anion2 --
    
    !---------------------------------- Somme Sur Les Especes Neutres --
    !--------------------------------------------------- using tLamda --
    tmp= Zero
    do jN=1,nSp
      if(vIonType(jN)/=0) cycle !jN is Neutral
      do J=1,nSp
        !
        X= Zero
        do K=J,nSp
          X= X +vMolal(K)*tZeta(J,K,jN)
        enddo
        !
        tmp= tmp + vMolal(jN)*vMolal(J) *(tLamda(J,jN) +X)
        !
      enddo
    end do
    if(iDebug>2) write(ff,'(I3,A,G15.6)') jN," fOsmo, Neutres=", tmp
    fOsmo= fOsmo + tmp
    !--------------------------------------------------------/neutres --
    !
    fOsmo= fOsmo *2.0D0
    if(iDebug>2) write(ff,'(I3,A,G15.6)') jN," fOsmo=", fOsmo
    !
    Osmo= One + 2.0D0*fOsmo/M
    
  else ! Force ionique=  0
  
    if (M==Zero) then
      Osmo=  One
    else
      Osmo=  One + fOsmo/M
    end if
  
  end if
  !
  !--------------------------------------------/COEFFICIENT OSMOTIQUE --
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  !------------------------------------------ COEFFICIENTS D'ACTIVITE --
  !
  tmp= 2.303D0 *0.5092D0 / 3.D0
  Fi1= -4.0d0*tmp      /1.2D0*log(One + 1.2D0*RI) &
  &    -2.0d0*tmp*RI/(One + 1.2D0*RI)
  BigF= 0.5d0 *Fi1
  ! BigF= BigF
  ! + SUM( m_c  m_a  B'_ca    ) 
  ! + SUM( m_c1 m_c2 PHI'_c1c2) 
  ! + SUM( m_a1 m_a2 PHI'_a1a2) 
  do I=1,nSp
    do J=1,nSp
      if(I==J) cycle
      tmp= vMolal(I)*vMolal(J)
      if(vIonType(I)*vIonType(J)== 1) BigF= BigF + Tmp *tETheta1(I,J)
      if(vIonType(I)*vIonType(J)==-1) BigF= BigF + Tmp *Bmx2(I,J)
    end do
  end do
  if(iDebug>2) write(ff,'(A,G15.6)') " BigF=", BigF
  
  !-------------------------------------------------- CALCUL LN GAMMA --
  do iC=1,nSp
    !
    !---------------------------------------- CALCUL POUR LES CATIONS --
    if(vIonType(iC)==1) then
      !
      AC= 0.5d0 *vCharge(iC)**2 *Fi1
      AC= vCharge(iC)**2 *BigF
      !
      if(iDebug>2) write(ff,'(/,A,/)') trim(vSpc(iC)%NamSp)
      if(iDebug>2) write(ff,'(A,G15.6)') " F-term=", AC
      !
      !----------------------------------------- Somme Sur Les Anions --
      Tmp= Zero
      do jA=1,nSp
        if(vIonType(jA)==-1) then !-> anion
          K1= MIN(iC,jA)
          K2= MAX(iC,jA)
          Tmp= Tmp + vMolal(jA) &
          & *( 2.0d0*Bmx1(K1,K2) +BigZ*Cmx(K1,K2) )
        end if
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 1=anions_=", Tmp
      AC= AC + Tmp
      !
      !---------------------------------------- Somme Sur Les Cations --
      Tmp= Zero
      do jC=1,nSp
        if(iC==jC) cycle
        if(vIonType(jC)==1) then !-> cation
          !
          K1= MIN(iC,jC)
          K2= MAX(iC,jC)
          !
          X= Zero
          do jA=1,nSp
            if(vIonType(jA)==-1) & ! iC-jC-jA
            & X= X &
            &  + vMolal(jA) *tPsi(K1,K2,jA)
          end do
          X= Zero
          !
          Tmp= Tmp &
          &  + vMolal(jC) &
          &    *( 2.0D0*tTheta(K1,K2) + 2.0D0*tETheta(K1,K2) + X)
          
          if(iDebug>2) write(FF,'(2A12,G12.3)') &
          & trim(vSpc(K1)%NamSp),trim(vSpc(K2)%NamSp), &
          & tTheta(K1,K2)+tETheta(K1,K2)
          
          !! Tmp= Tmp &
          !! &  + vMolal(jC) &
          !! &    *( 2.0D0*tTheta(jC,iC) + 2.0D0*tETheta(jC,iC) )
          !
        end if
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 2=cations=", Tmp
      AC= AC + Tmp
      !
      !------------------------------------------ Somme Anion - Anion --
      Tmp= Zero
      do jA=1,nSp-1
        if(vIonType(jA)==-1) then !ANI
          do kA=jA+1,nSp
            if(vIonType(jA)==-1) then !ANI
              Tmp= Tmp &
              & + vMolal(jA) *vMolal(kA) *tPsi(jA,kA,iC) !ANI-ANI-CAT
            end if
          enddo
        end if
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 3=ani/ani=", Tmp
      AC=  AC + Tmp
      !
      !------------------ Somme Croisée Sur Les Cations Et Les Anions --
      Tmp=  Zero
      do I=1,nSp
        do J=1,nSp
          if(vIonType(I)*vIonType(J)==-1) &
          & Tmp= Tmp &
          &    + vMolal(I)*vMolal(J) *vCharge(iC) *Cmx(I,J)
        enddo
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 4=ani/cat=", Tmp
      AC= AC + Tmp
      !
      !-------------------------------- Somme Sur Les Especes Neutres --
      Tmp= Zero
      do jN=1,nSp
        if(vIonType(jN)==0) &
        & Tmp= Tmp +2.0d0*vMolal(jN)*tLamda(iC,jN)
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 5=neutral=", Tmp
      AC= AC +Tmp
      !---/
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      vLnGam(iC)=AC
    end if
    !---------------------------------------/ CALCUL POUR LES CATIONS --
  end do
    
  do iA=1,nSp
    !
    !----------------------------------------- CALCUL POUR LES ANIONS --
    if(vIonType(iA)==-1) then
      !
      if(iDebug>2) write(ff,'(/,A,/)') trim(vSpc(iA)%NamSp)
      !
      AC= 0.5d0 *vCharge(iA)**2 *Fi1
      AC= vCharge(iA)**2 *BigF
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      !-------------------------------------------- somme sur cations --
      Tmp= Zero
      do jC=1,nSp
        if(vIonType(jC)==1) then
          K1= MIN(iA,jC)
          K2= MAX(iA,jC)
          Tmp= Tmp + vMolal(jC) &
          &         *(2.0d0 *Bmx1(K1,K2) &
          &            + BigZ*Cmx(K1,K2))
        end if
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 1=cations=", Tmp
      AC= AC + Tmp
      !
      !--------------------------------------------- somme sur anions --
      Tmp= Zero
      do jA=1,nSp
        if(iA==jA) cycle
        if(vIonType(jA)==-1) then
          !
          K1= MIN(iA,jA)
          K2= MAX(iA,jA)
          !
          X= Zero
          do jC=1,nSp
            if(vIonType(jC)==1) &
            & X= X &
            &  + vMolal(jC) *tPsi(K1,K2,jC)
          end do
          !X= Zero
          Tmp=  Tmp &
          &  + 2.0d0*vMolal(jA)*(tTheta(K1,K2) +tETheta(K1,K2)) 
          !
          if(iDebug>2) write(FF,'(2A12,G12.3)') &
          & trim(vSpc(K1)%NamSp),trim(vSpc(K2)%NamSp), &
          & tTheta(K1,K2)+tETheta(K1,K2)
          !
        endIf
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 2=anions =", Tmp
      AC= AC + Tmp
      !
      !-------------------------------------- somme cations - cations --
      Tmp= Zero
      do I=1,nSp-1
        if(vIonType(I)==1) then !CAT
          do J=I+1,nSp
            if(vIonType(J)==1) then !CAT
              !
              Tmp= Tmp &
              &  + vMolal(I)*vMolal(J) *tPsi(I,J,iA) !CAT-CAT-ANI
              !
            end if
          enddo
        end if
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 3=cat/cat=", Tmp
      AC= AC + Tmp
      !
      !------------------ somme croisée sur les cations et les anions --
      Tmp= Zero
      do I=1,nSp
        do J=1,nSp
          if(vIonType(I)*vIonType(J)==-1) &
          & Tmp= Tmp &
          &    - vMolal(I)*vMolal(J) *vCharge(iA) *Cmx (I,J)
        enddo
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " 4=ani/cat=", Tmp
      AC=  AC + Tmp
      !
      !------------------------------------------ Somme Anion - Anion --
      !! Tmp= Zero
      !! do jA=1,nSp-1
      !!   if(vCharge(jA)<0) then
      !!     !! do kA=jA+2,nSp
      !!     do kA=jA+1,nSp
      !!       if(vCharge(kA)<0) &
      !!       & Tmp= Tmp + vMolal(jA)*vMolal(kA)*tETheta1(jA,kA)
      !!     enddo
      !!   end if
      !! enddo
      !! Tmp= Tmp *vCharge(iA)**2
      !! if(iDebug>2) write(ff,'(A,G15.6)') " ani/ani=", Tmp
      !! AC= AC + Tmp
      
      !-------------------------------- Somme Sur Les Especes Neutres --
      Tmp= Zero
      do jN=1,nSp
        if(vIonType(jN)==0) &
        & Tmp= Tmp +2.0d0*vMolal(jN)*tLamda(iA,jN)
      enddo
      !! Tmp= SUM(vMolal(:)*tLamda(iA,:),MASK=vIonType(:)==0)
      if(iDebug>2) write(ff,'(A,G15.6)') " 5=neutral=", Tmp
      AC= AC + Tmp
      !---/
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      vLnGam(iA)=AC
    end if
    !----------------------------------------/ CALCUL POUR LES ANIONS --
  end do
  
  do jN=1,nSp
    !
    !-------------------------------- CALCUL POUR LES ESPECES NEUTRES --
    if(vIonType(jN)==0) then
      !
      if(iDebug>2) write(ff,'(/,A,/)') trim(vSpc(jN)%NamSp)
      !
      AC= Zero
      !
      Tmp= Zero
      do J=1,nSp
        Tmp= Tmp + 2.0*vMolal(J)*tLamda(J,jN)
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " tLamda_=", Tmp
      AC= AC + Tmp
      !
      !----------------------------------------------- terme ternaire --
      Tmp= Zero
      do J=1,nSp
        if(vCharge(J)/=0) then
          do K=J+1,nSp
            if(vCharge(K)/=0) then
              Tmp= Tmp + vMolal(J)*vMolal(K)*tZeta(J,K,jN)
            end if
          enddo
        end if
      enddo
      if(iDebug>2) write(ff,'(A,G15.6)') " ternary=", Tmp
      !-----------------------------------------------/terme ternaire --
      !
      AC= AC + Tmp
      !
      vLnGam(jN)=AC
    end if
    !-------------------------------/ CALCUL POUR LES ESPECES NEUTRES --
    !
  enddo
  !
  do I=1,nSp
    if(iDebug>2) &
    & write(ff,'(I3,1X,A15,G15.6)') I,vSpc(I)%NamSp,exp(vLnGam(I))
  enddo
  !
  if(iDebug>2) close(ff)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Pitzer_Calc_Gamma"
  
  return

contains

subroutine CheckCoeffs

  write(FF,'(A)') "==Beta0-1-2 Cphi=="
  do I=1,nSp
    do J=1,nSp
      if(J>I) &
      & write(FF,'(2A12,4G12.3)') &
      & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp), &
      & tBeta0(I,J), tBeta1(I,J), tBeta2(I,J), tCPhi(I,J)
    end do
  end do
  
  write(FF,'(A)') "==Psi, CAT-CAT-ANI=="
  do K=1,nSp
    do I=1,nSp-1
      do J=I+1,nSp
        if( vIonType(I)== 1 .and. &
        &   vIonType(J)== 1 .and. &
        &   vIonType(K)==-1 )     &
        & write(FF,'(3A12,G12.3)') &
        & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp),trim(vSpc(K)%NamSp), &
        & tPsi(I,J,K)
      end do
    end do
  end do
  
  write(FF,'(A)') "==Psi, ANI-ANI-CAT=="
  do K=1,nSp
    do I=1,nSp-1
      do J=I+1,nSp
        if( vIonType(I)==-1 .and. &
        &   vIonType(J)==-1 .and. &
        &   vIonType(K)== 1 )     &
        & write(FF,'(3A12,G12.3)') &
        & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp),trim(vSpc(K)%NamSp), &
        & tPsi(I,J,K)
      end do
    end do
  end do
  
  write(FF,'(A)') "==Lambda, ION-NEU=="
  do J=1,nSp
    if( vIonType(J)/=0 ) cycle
    do I=1,nSp
      if( vIonType(I)/=0 )     &
      & write(FF,'(2A12,G12.3)') &
      & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp), &
      & tLamda(I,J)
    end do
  end do
  
  !! write(FF,'(A)') "==Beta0-1-2 Cphi=="
  !! do I=1,nSp
  !!   do J=1,nSp
  !!     if(J<I) &
  !!     & write(FF,'(2A12,4G12.3)') &
  !!     & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp), &
  !!     & tBeta0(I,J), tBeta1(I,J), tBeta2(I,J), tCPhi(I,J)
  !!   end do
  !! end do
  !! 
  !! write(FF,'(A)') "==Psi=="
  !! do K=1,nSp
  !!   do I=1,nSp
  !!     do J=1,nSp
  !!       if( J<I &
  !!       .and. vIonType(I)== 1  &
  !!       .and. vIonType(J)== 1  &
  !!       .and. vIonType(K)==-1) &
  !!       & write(FF,'(3A12,G12.3)') &
  !!       & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp),trim(vSpc(K)%NamSp), &
  !!       & tPsi(I,J,K)
  !!     end do
  !!   end do
  !! end do
  !! do K=1,nSp
  !!   do I=1,nSp
  !!     do J=1,nSp
  !!       if( J<I &
  !!       .and. vIonType(I)==-1  &
  !!       .and. vIonType(J)==-1  &
  !!       .and. vIonType(K)== 1) &
  !!       & write(FF,'(3A12,G12.3)') &
  !!       & trim(vSpc(I)%NamSp),trim(vSpc(J)%NamSp),trim(vSpc(K)%NamSp), &
  !!       & tPsi(I,J,K)
  !!     end do
  !!   end do
  !! end do
  
  return
end subroutine CheckCoeffs

end subroutine Pitzer_Calc_Gamma

real(dp) function FONC1(X)
!--
!-- function g, used for computing B_MX and B_MX'
!-- from Beta0, Beta1, Beta2
!--
  real(dp)::X
  if (X==Zero) then  ;  FONC1= Zero
  else               ;  FONC1= 2.0d0*(One - (One + X)*exp(-X)) /X/X
  end if
end function FONC1

real(dp) function FONC2(X)
!--
!-- function g', used for computing B_MX and B_MX'
!-- from Beta0, Beta1, Beta2
!--
  real(dp)::X
  if (X==Zero) then  ;  FONC2= Zero
  else               ;  FONC2=-2.0d0*(One - (One + X + X*X/2.0d0)*exp(-X)) /X/X
  end if
end function FONC2

subroutine EThetaCalc ( &
& Zi, Zj, SqrtI, Aphi,  &
& Eth, Eth1)
!--
!-- Calcul de la fonction ETheta(I) --
!-- (effet electrostatique asymetrique) --
!-- par la theorie de Pitzer (1975) --
!-- IN --
!  Zi: charge du premier ion
!  Zj: charge du deuxieme ion
!  SqrtI: racine carree de la force ionique
!  Aphi:  terme de Debye-Huckel
!-- OUT --
!  Eth : effet electrostatique asymetrique
!  Eth1 : derivee de Eth par rapport a la temperature
!--
  integer, intent(in) :: ZI, ZJ
  real(dp),intent(in) :: SqrtI, APHI
  real(dp),intent(out):: eth, eth1
  !
  real(dp):: J0mn,J1mn,J0mm,J1mm,J0nn,J1nn
  real(dp):: XMN,XMM,XNN,RI2
  !
  eth=  Zero
  eth1= Zero
  RI2=  SqrtI**2
  
  if (SqrtI==Zero .or. Zi==Zj) return
  
  Xmn= 6.0d0 *Zi*Zj *Aphi *SqrtI
  Xmm= 6.0d0 *Zi*Zi *Aphi *SqrtI
  Xnn= 6.0d0 *Zj*Zj *Aphi *SqrtI
  
  call EThetaParam (Xmn, J0mn, J1mn)
  call EThetaParam (Xmm, J0mm, J1mm)
  call EThetaParam (Xnn, J0nn, J1nn)
  
  Eth=  Zi*Zj *(J0mn     - 0.5d0*(    J0mm +     J0nn)) /4.D0/RI2
  Eth1= Zi*Zj *(Xmn*J1mn - 0.5d0*(Xmm*J1mm + Xnn*J1nn)) /8.D0/RI2/RI2 &
  &   - eth/RI2
  
  return
end subroutine EThetaCalc

subroutine EThetaParam (X, JX, J1X)
!--
!-- Element pour le calcul de la fonction ETheta(I)
!-- (effet electrostatique asymetrique)
!-- par la theorie de Pitzer (1975)
!-- IN-
!--   X : charge de l'ion considere
!-- OUT-
!--   JX : intermediaire de calcul
!--   J1X : intermediaire de calcul

  real(dp),intent(in) :: X
  real(dp),intent(out):: JX, J1X
  
  real(dp):: Z, DZ
  integer :: I
  real(dp),dimension(0:22):: aik, aiik, bk, dk
  !
  data (aik(i), aiik(i), i=0,22)                &
  & /1.925154014814667D0,  0.628023320520852D0, &
  & -0.060076477753119D0,  0.462762985338493D0, &
  & -0.029779077456514D0,  0.150044637187895D0, &
  & -0.007299499690937D0, -0.028796057604906D0, &
  &  0.000388260636404D0, -0.036552745910311D0, &
  &  0.000636874599598D0, -0.001668087945272D0, &
  &  0.000036583601823D0,  0.006519840398744D0, &
  & -0.000045036975204D0,  0.001130378079086D0, &
  & -0.000004537895710D0, -0.000887171310131D0, &
  &  0.000002937706971D0, -0.000242107641309D0, &
  &  0.000000396566462D0,  0.000087294451594D0, &
  & -0.000000202099617D0,  0.000034682122751D0, &
  & -0.000000025267769D0, -0.000004583768938D0, &
  &  0.000000013522610D0, -0.000003548684306D0, &
  &  0.000000001229405D0, -0.000000250453880D0, &
  & -0.000000000821969D0,  0.000000216991779D0, &
  & -0.000000000050847D0,  0.000000080779570D0, &
  &  0.000000000046333D0,  0.000000004558555D0, &
  &  0.000000000001943D0, -0.000000006944757D0, &
  & -0.000000000002563D0, -0.000000002849257D0, &
  & -0.000000000010991D0,  0.000000000237816D0, &
  &  Zero,                Zero,               &
  &  Zero,                Zero                /
  
  do i=20,0,-1
  
    if (x <= One) then
      Z=     4.0D0 * x**0.2D0    - 2.0D0
      dZ=    0.8D0 * x**(-0.8D0)
      bk(i)= Z*bk(i+1) - bk(i+2) + aik(i)
    else
      Z=     40.0D0/9.0D0 *x**(-0.1D0) - 22.0D0/9.0D0
      dZ=   -40.0D0/90.0D0*x**(-1.1D0)
      bk(i)= Z*bk(i+1) - bk(i+2) + aiik(i)
    end if
    dk(i)=  bk(i+1) + Z*dk(i+1) - dk(i+2)
  
  enddo
  
  JX=  0.25D0*x - One + 0.5D0*(bk(0) - bk(2))
  J1X= 0.25D0   +    0.5D0*dZ*(dk(0) - dk(2))
  
  return
end subroutine EThetaParam

subroutine Pitzer_Calc_APhi(RhoW,EpsW,TdgK,APhi)
!--
!-- Aphi calculation is slightly different from that used in DH
!-- refs: Anderson-Crerar, p449; or Lin_2003_FPE ...
!--
  ! use M_Numeric_Const,only: Pi
  ! use M_Dtb_Const,only: Electron,Boltzmann,Avogadro,EpsVacuum
  !
  real(dp),intent(in) :: RhoW,EpsW,TdgK
  real(dp),intent(out):: APhi
  
  !! real(dp):: X
  
  !! X= SQRT( (2.0D0*pi*Avogadro *1.0D-3)  &
  !! &          /(4.0D0*pi *EpsVacuum *Boltzmann)**3 ) &
  !! &   * Electron**3 &
  !! &   *1.0D+3 /3.0D0
  !! print '(F24.6)',X
  !! pause
  
  !! Aphi= SQRT( (2.0D0*pi*Avogadro *RhoW*1.0D-3)  &
  !! &          /(4.0D0*pi *EpsVacuum*EpsW *Boltzmann*TdgK)**3 ) &
  !! &   * Electron**3 &
  !! &   *1.0D+3 /3.0D0
  
  Aphi= 1400684.314D0 *SQRT( RhoW / (EpsW*TdgK)**3 )
  
  !~ APhi_DH= 1.824829238D6 *SQRT( RhoW / (EpsW*TdgK)**3)
  !~ call Pitzer_Calc_APhiMonnin(TdgK,1.0D0,Aphi_)
  !~ print '(A,3G15.6)',"Aphi, AphiMonnin,APhi_DH:", Aphi,Aphi_,APhi_DH
  !~ pause
  
  return
end subroutine Pitzer_Calc_APhi

subroutine Pitzer_Calc_APhiMonnin(TdgK,Pbar,Aphi_)
!--
!-- Module de calcul du terme APhi de Debye-Huckel
!-- (version spéciale Pitzer)
!-- Correlations de Monnin, Chemical Geology, 153(1999) 187-209
!--
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(out):: Aphi_
  !
  real(dp)::Psat,T,A0,AV
  !
  A0=Zero
  AV=Zero
  APhi_=Zero
  T=TdgK
  !
  A0=  3.36901531D-1                 &
  &  - 6.32100430D-4 * T             &
  &  + 9.14252359    / T             & 
  &  - 1.35143986D-2 * log(T)        & 
  &  + 2.26089488D-3 / (T-263.0D0)   &
  &  + 1.92118597D-6 * T*T           &
  &  + 4.52586464D1  / (680.D0-T)     
  !
  if (T<=373.15D0) then
    !
    AV=  8.106377         &
    & - 1.256008D-1 *T    &
    & + 7.760276D-4 *T**2 &
    & - 2.098163D-6 *T**3 &
    & + 2.25777D-9  *T**4 !cm3.kg^0.5.mol(-1.5), Monnin,1989
    !
    AV=  AV * 1.0D-06 !conversion en m3
    !
  else
    !
    AV= 3.849971D2          &
    & - 6.982754     * T    &
    & + 3.877068D-2  * T**2 &
    & - 1.11381D-4   * T**3 &
    & + 1.589736D-7  * T**4 &
    & - 9.395266D-11 * T**5 &
    & + 6.469241D4 / (680.D0-T) !cm3.kg^0.5.mol(-1.5), Monnin,1989
    !
    AV=  AV * 1.0D-06 !conversion en m3
    !
  end if
  !  
  Psat= exp( 7.3649D1 -7.2582D3/T -7.3037D0*log(T) +4.1653D-6*T*T )
  ! in Pascal
  !
  Aphi_ =  A0 + (Pbar*1.0D+5 - Psat) * (-AV/(4.0*8.314*T))  !Av=-1/4RT * dAphi/dP_, Monnin, 1989
  !  
  return
end subroutine Pitzer_Calc_APhiMonnin

end module M_Solmodel_Calc_Pitzer
