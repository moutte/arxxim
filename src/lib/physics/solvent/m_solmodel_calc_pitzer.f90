module M_Solmodel_Calc_Pitzer
!--
!-- Module routines associées à Pitzer
!-- Livraison 10/2006: Nicolas Ferrando
!-- Remoduling : Anthony Michel
!-- MODIFIED : J.Moutte, 11/2006
!-- -> M_Solmodel_Pitzer_Dtb=  database reading, build vPitz, vIPitz
!-- -> M_Solmodel_Pitzer_Calc= calculations
!--

  use M_Kinds
  use M_Trace,only: fTrc,iDebug,T_
  
  implicit none
  
  private
  !
  public:: Solmodel_Calc_Pitzer
  public:: EThetaCalc
  
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
  type(T_Species),dimension(:),allocatable:: vSolut
  real(dp),       dimension(:),allocatable:: vMolalSolut
  integer :: I,J,N
  real(dp):: APhi
  
  call Pitzer_Calc_APhi( &
  & Rho,Eps,TdgK, &
  & APhi)
  
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
  end do
  
  call Pitzer_Calc_Gamma( &
  & vSolut,Aphi,vMolalSolut, & !IN
  & Osmo,vLnGam) !OUT
  !
  LnActSv= -Osmo *MWsv * (sum(vMolalSolut))
  
  deallocate(vSolut)
  deallocate(vMolalSolut)
  
  return
end subroutine Solmodel_Calc_Pitzer

subroutine Pitzer_Calc_Gamma( &
& vSpc,Aphi,vMolal, &
& Osmo,vLnGam) !,ERR)
!--
!-- Calcule les coefficients d'activité des solute's
!-- et le coefficient osmotique de l'eau par le modele de Pitzer
!--
!-- vSpc is a list of aqu'species with (kanarazu) solvent as vSpc(1)
!-- and solutes as vSpc(2:nSp)
!--
  use M_T_Species,only: T_Species
  
  !IN
  !  vSpc:   solute's
  !  Aphi:   pente du terme de Debye-Huckel
  !  vMolal: molalite's des solute's
  !
  !  !vCharge: charge des solute's 
  !  !tBeta0: parametre du modele de Pitzer
  !  !tBeta1: parametre du modele de Pitzer.
  !  !tBeta2: parametre du modele de Pitzer
  !  !tCPhi:  parametre du modele de Pitzer
  !  !tTheta: parametre du modele de Pitzer
  !  !tPsi:   parametre du modele de Pitzer
  !  !tAlfa1: parametre du modele de Pitzer
  !  !tAlfa2: parametre du modele de Pitzer
  !OUT
  !  Osmo: Coefficient osmotique 
  !  vLnGam: log neperien des coefficient d'activite des constituants
  !  ERR: indicateur d'erreur (non utilise -> a completer si besoin)
  
  use M_Solmodel_Pitzer_Dtb,only: tBeta0,tBeta1,tBeta2,tAlfa1,tAlfa2
  use M_Solmodel_Pitzer_Dtb,only: tTheta,tPsi,tCPhi,tZeta,tLamda
  use M_IoTools,only:GetUnit
  !
  type(T_Species),dimension(:),intent(in) :: vSpc
  real(dp),                    intent(in) :: Aphi
  real(dp),       dimension(:),intent(in) :: vMolal
  real(dp),                    intent(out):: Osmo
  real(dp),       dimension(:),intent(out):: vLnGam
  !
  !--------------------------------------------------- Variables locales
  integer :: I,J,K,iC,iA,jC,jA,jN,kC,kA,nSp,A1,A2,C1,C2
  integer :: zA !, zC
  integer :: FF !file index
  real(dp):: Ionic, M, BigZ
  real(dp):: RI, Fi0, Fi1, fOsmo
  real(dp):: AC, Tmp
  real(dp):: ETheta, ETheta1
  real(dp),dimension(1:size(vSpc),1:size(vSpc)):: &
  & Bmx1, & !B_MX  of Pitzer equations, p104
  & Bmx2, & !B'_MX of Pitzer equations, p104
  & Cmx     !C_MX  of Pitzer equations, p104
  integer,allocatable:: vCharge(:),vIonType(:)
  !--------------------------------------------------------------------/
  !
  !Fonctions utilitaires
  !real(dp) FONC1
  !real(dp) FONC2
  !integer  IonType
  !
  if(idebug>1) write(fTrc,'(A)') "< Pitzer_Calc_Gamma"
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
    call GetUnit(ff)
    open(ff,file="debug_pitzer.log")
    write(ff,'(A,/)') "trace of routine Pitzer_Calc_Gamma"
  end if
  !ERR=0
  !
  !--------------------------------------------------------FORCE IONIQUE
  M=    Zero  !! total molality
  BigZ= Zero  !! total molal charge, sum(|z_i|.m_i)
  Ionic=Zero  !! ionic strength
  do I=1,nSp
    M=     M     + vMolal(I)
    BigZ=  BigZ  + vMolal(I)*vCharge(I)*IonType(I,1) !-> sum(|z_i|.m_i)= BigZ
    Ionic= Ionic + vMolal(I)*vCharge(I)*vCharge(I)
  end do
  Ionic= Ionic/2.0D0
  if(iDebug>2) write(ff,'(A,G15.6)') "Ionic= ", Ionic
  !-------------------------------------------------------/FORCE IONIQUE
  !
  !-------------------------------------------COEFFICIENTS LONGUE PORTEE
  RI=  SQRT(Ionic)
  Fi0= -4.0d0*Aphi*Ionic/1.2D0*log(One + 1.2D0*RI)
  Fi1= -4.0d0*Aphi      /1.2D0*log(One + 1.2D0*RI) &
  &    -2.0d0*Aphi*RI/(One + 1.2D0*RI)
  if(iDebug>2) write(ff,'(A,G15.6)') "Fi0=   ", Fi0
  if(iDebug>2) write(ff,'(A,G15.6)') "Fi1=   ", Fi1
  !------------------------------------------/COEFFICIENTS LONGUE PORTEE
  !
  Bmx1= Zero
  Bmx2= Zero
  Cmx=  Zero
  !
  if (Ionic /= Zero) then
    !----------------------------------------PARAMETRES DE PITZER B et C
    !--------------------------------using tBeta0, tBeta1, tBeta2, tCPhi
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
        if (vCharge(I)*vCharge(J)/=0) then !C_cation-anion eq24, Moeller'98
          Cmx(I,J)=  0.5d0*tCPhi(I,J) / (abs(vCharge(I)*vCharge(J)))**0.5
        else
          Cmx(I,J)=  Zero
        end if
        !
        if(iDebug>2 .and. Bmx1(I,J)/=Zero) write(ff,'(A,3(G15.6,1X),2(A,A1))') &
        & "Bmx1,Bmx2,Cmx", &
        & Bmx1(I,J), Bmx2(I,J), Cmx(I,J), &
        & trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_
        !
      end do
    end do
    !---------------------------------------/PARAMETRES DE PITZER B et C
    !
    !----------------------------------------------COEFFICIENT OSMOTIQUE
    fOsmo=  Ionic*Fi1 - Fi0
    !
    !-----------------------------------------somme croisee anion/cation
    do I=1,nSp
      do J=1,nSp
        if(IonType(I,1)*IonType(J,-1)/=0) then !J(+)I(-) or J(-)I(+)
          fOsmo= fOsmo &
          &    + 2.0d0*vMolal(I)*vMolal(J) &
          &      * ( Bmx1(I,J) +Ionic*Bmx2(I,J) + 2.0d0*BigZ*Cmx(I,J) )
        end if
      end do
      !
      if(iDebug>2) write(ff,'(I3,A,G15.6)') I," fOsmo, anion/cation=", fOsmo
      !call Pause_
    end do
    !--------------------------------------somme croisee Cation1/Cation2
    
    !-------------------------------------------------using tPsi, tTheta
    do C1=1,nSp
      do C2=C1,nSp
        if(IonType(C1,1)*IonType(C2,1)*(C1-C2)/=0) then
          Tmp=  Zero
          do K=1,nSp !if K is anion
            !N.Ferrando considers K is anion or neutral
            !Tmp= Tmp + vMolal(K)*tPsi(C1,C2,K) *(max(IonType(K,-1),IonType(K,0)))
            Tmp= Tmp + vMolal(K) *tPsi(C1,C2,K) *IonType(K,-1)
          end do
          call EThetaCalc(vCharge(C1),vCharge(C2), RI, Aphi, ETheta, ETheta1)
          fOsmo= fOsmo &
          &    + vMolal(C1)*vMolal(C2) &
          &      *(tTheta(C1,C2) +ETheta +RI*RI*ETheta1 +Tmp)
        end if
      end do
      !
      if(iDebug>2) write(ff,'(I3,A,G15.6)') C1," fOsmo, cation1/cation2=", fOsmo
      !call Pause_
    end do
    
    !--=================================< Somme Croisee Anion1/Anion2 ==
    !--==========================================< using tPsi, tTheta ==
    do A1=1,nSp
      do A2=A1,nSp
        if(IonType(A1,-1)*IonType(A2,-1)*(A1-A2)/=0) then
          !A1/=A2 & A1(-) & A2(-)
          Tmp=  Zero
          do K=1,nSp !if K is cation
            !N.Ferrando considers K is cation or neutral
            !Tmp= Tmp + vMolal(K) *tPsi(A1,A2,K) *MAX(IonType(K,1),IonType(K,0))
            Tmp= Tmp + vMolal(K) *tPsi(A1,A2,K) *IonType(K,1)
          end do
          call EThetaCalc(vCharge(A1),vCharge(A2), RI, Aphi, ETheta, ETheta1)
          fOsmo= fOsmo &
          &    + vMolal(A1) *vMolal(A2) *(tTheta(A1,A2) +ETheta +RI*RI*ETheta1 +Tmp)
        end if
      end do
      if(iDebug>2) write(ff,'(I3,A,G15.6)') A1," fOsmo, Anion1/Anion2=", fOsmo
      !call Pause_
    end do
    !
    !--===============================< Somme Sur Les Especes Neutres ==
    !==================================================< using tLamda ==
    do jN=1,nSp
      if(IonType(jN,0)/=0) then !jN Is Neutral
        do J=1,nSp
          Tmp=Zero
          do K=J,nSp
            Tmp= Tmp +vMolal(K)*tZeta(J,K,jN)
          end do
          fOsmo= fOsmo + vMolal(jN)*vMolal(J) *(tLamda(J,jN) +Tmp)
        end do
      end if
      if(iDebug>2) write(ff,'(I3,A,G15.6)') jN," fOsmo, Neutres=", fOsmo
      !call Pause_
    end do
    !
    Osmo= One + fOsmo/M
    !
  else ! Force ionique=  0
  
    if (M==Zero) then
      Osmo=  One
    else
      Osmo=  One + fOsmo/M
    end if
  
  end if
  !
  !--=======================================< COEFFICIENTS D'ACTIVITE ==
  !
  !--=======================================< CALCUL POUR LES CATIONS ==
  do iC=1,nSp
    if(vCharge(iC)>0) then
      !
      if(iDebug>2) write(ff,'(/,A,/)') trim(vSpc(iC)%NamSp)
      !
      AC= 0.5d0 *vCharge(iC)**2 *Fi1
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      !--======================================< Somme Sur Les Anions ==
      Tmp= Zero
      do jA=1,nSp
        if(vCharge(jA)<0) &
        !Tmp=  Tmp + vMolal(jA) *(Bmx1(iC,jA) +BigZ*Cmx(iC,jA)) *IonType(jA,-1)
        Tmp= Tmp + 2.0d0*vMolal(jA) *( Bmx1(iC,jA) +BigZ*Cmx(iC,jA) )
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " anions_=", Tmp
      AC= AC + Tmp
      !
      !--=====================================< Somme Sur Les Cations ==
      Tmp=  Zero
      do jC=1,nSp
        if(vCharge(jC)>0) then
          call EThetaCalc(vCharge(iC),vCharge(jC), RI, Aphi, ETheta, ETheta1)
          !Tmp= Tmp + vMolal(jC)*( tTheta(min(iC,jC),MAX(iC,jC)) + ETheta )
          Tmp= Tmp + 2.0d0*vMolal(jC)*( tTheta(min(iC,jC),MAX(iC,jC)) + ETheta )
        end if
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " cations=", Tmp
      AC= AC + Tmp
      !
      !--===============< Somme Croisée Sur Les Cations Et Les Anions ==
      Tmp=  Zero
      do J=1,nSp
        do K=1,nSp !when K(-) and J(+)
          Tmp= Tmp &
          &  + vMolal(J)*vMolal(K) &
          &    * ( vCharge(iC)**2 *Bmx2(J,K)   &
          &      + vCharge(iC)    *Cmx(J,K)    &
          &      + tPsi(min(iC,J),MAX(iC,J),K) &
          &      ) &
          &    * IonType(J,1)*IonType(K,-1)
        end do
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " ani/cat=", Tmp
      AC= AC + Tmp
      !
      !--=======================================< Somme Anion < Anion ==
      Tmp= Zero
      do jA=1,nSp-1
        if(vCharge(jA)<0) then
          do kA=jA+1,nSp
            if(vCharge(kA)<0) then
              call EThetaCalc(vCharge(jA),vCharge(kA), RI, Aphi, ETheta, ETheta1)
              Tmp= Tmp &
              & + vMolal(jA) *vMolal(kA) *( vCharge(iC)**2 *ETheta1 + tPsi(jA,kA,iC) )
            end if
          end do
        end if
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " ani/ani=", Tmp
      AC=  AC + Tmp
      !
      !--< Somme Cations < Cations
      Tmp= Zero
      do jC=1,nSp-1
        if(vCharge(jC)>0) then
          !do kC=jC+2,nSp
          do kC=jC+1,nSp
            if(vCharge(kC)>0) then
              call EThetaCalc(vCharge(jC),vCharge(kC), RI, Aphi, ETheta, ETheta1)
              Tmp= Tmp + vMolal(jC)*vMolal(kC) *ETheta1
            end if
          end do
        end if
      end do
      Tmp= Tmp *vCharge(iC)**2
      if(iDebug>2) write(ff,'(A,G15.6)') " cat/cat=", Tmp
      AC=  AC + Tmp
      !
      !--< Somme Sur Les Especes Neutres
      Tmp= Zero
      do jN=1,nSp
        Tmp= Tmp +2.0d0*vMolal(jN)*tLamda(iC,jN)*IonType(jN,0)
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " neutral=", Tmp
      AC= AC +Tmp
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      vLnGam(iC)=AC
    end if
  end do
  !--======================================</ CALCUL POUR LES CATIONS ==
  
  !--========================================< CALCUL POUR LES ANIONS ==
  do iA=1,nSp
    
    if(vCharge(iA)<0) then
      !
      if(iDebug>2) write(ff,'(/,A,/)') trim(vSpc(iA)%NamSp)
      !
      zA= vCharge(iA)
      !
      AC= 0.5d0 *zA**2 *Fi1  !*IonType(iA,-1)
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      !--=========================================< somme sur cations ==
      Tmp= Zero
      do jC=1,nSp
        Tmp= Tmp + 2.0d0 *vMolal(jC) *(Bmx1(jC,iA) + BigZ*Cmx(jC,iA)) *IonType(jC,1)
      end do
      AC= AC + Tmp
      if(iDebug>2) write(ff,'(A,G15.6)') " cations=", Tmp
      !
      !--==========================================< somme sur anions ==
      Tmp= Zero
      do jA=1,nSp
        !ZI=  vCharge(iA)
        !ZJ=  vCharge(jA)
        if(vCharge(jA)<0) then
          call EThetaCalc(zA,vCharge(jA), RI, Aphi, ETheta, ETheta1)
          Tmp=  Tmp + 2.0d0*vMolal(jA)*(tTheta(min(iA,jA),MAX(iA,jA)) +ETheta) 
        endIf
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " anions =", Tmp
      AC= AC + Tmp
      !
      !--===============< somme croisée sur les cations et les anions ==
      Tmp= Zero
      do J=1,nSp
        do K=1,nSp
          Tmp= Tmp &
          &  + vMolal(J)*vMolal(K) &
          &    * ( zA**2   *Bmx2(J,K)   &
          &      + abs(zA) *Cmx (J,K)    &
          &      + tPsi(min(iA,K),MAX(iA,K),J)   &
          &      ) &
          &    *IonType(J,1) *IonType(K,-1)
        end do
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " ani/cat=", Tmp
      AC=  AC + Tmp
      !
      !--===================================< somme cations < cations ==
      Tmp= Zero
      do jC=1,nSp-1
        if(vCharge(jC)>0) then
          do kC=jC+1,nSp !§§§!kC=jC+2,nSp 
            if(vCharge(kC)>0) then
              call EThetaCalc(vCharge(jC),vCharge(kC), RI, Aphi, ETheta, ETheta1)
              Tmp= Tmp &
              &  + vMolal(jC)*vMolal(kC) *( zA**2 *ETheta1 + tPsi(jC,kC,iA) )
            end if
          end do
        end if
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " cat/cat=", Tmp
      AC= AC + Tmp
      !
      !--=======================================< Somme Anion < Anion ==
      Tmp= Zero
      do jA=1,nSp-1
        if(vCharge(jA)<0) then
          !! do kA=jA+2,nSp
          do kA=jA+1,nSp
            !ZJ=  vCharge(jA)
            !ZK=  vCharge(kA)
            if(vCharge(kA)<0) then
              call EThetaCalc(vCharge(jA),vCharge(kA), RI, Aphi, ETheta, ETheta1)
              Tmp= Tmp + vMolal(jA)*vMolal(kA)*ETheta1
            end if
          end do
        end if
      end do
      Tmp= Tmp *vCharge(iA)**2
      if(iDebug>2) write(ff,'(A,G15.6)') " ani/ani=", Tmp
      AC= AC + Tmp
      
      !--=============================< Somme Sur Les Especes Neutres ==
      Tmp= Zero
      do jN=1,nSp
        !Tmp= Tmp + vMolal(jN)*Bmx1(iC,jN)*IonType(jN,0)
        Tmp= Tmp +2.0d0*vMolal(jN)*tLamda(iA,jN)*IonType(jN,0)
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " neutral=", Tmp
      AC= AC + Tmp
      !Tmp=  Zero
      !do J=1,nSp
      !   TMP1=Bmx1(iA,J)       
      !   Tmp= Tmp + vMolal(J)*TMP1*IonType(J,0)
      !end do
      !AC=  AC + 2.0d0*Tmp
      !
      if(iDebug>2) write(ff,'(A,G15.6)') " LnGamma=", AC
      !
      vLnGam(iA)=AC
    end if
    
  end do
  !--=======================================</ CALCUL POUR LES ANIONS ==
  !
  !--===============================< CALCUL POUR LES ESPECES NEUTRES ==
  do jN=1,nSp
  
    if(vCharge(jN)==0) then
      !
      if(iDebug>2) write(ff,'(/,A,/)') trim(vSpc(jN)%NamSp)
      !
      AC= Zero
      !
      Tmp= Zero
      do J=1,nSp
        Tmp= Tmp + 2.0*vMolal(J)*tLamda(J,jN)
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " tLamda_=", Tmp
      AC= AC + Tmp
      !
      !--< terme ternaire
      Tmp= Zero
      do J=1,nSp
        if(vCharge(J)/=0) then
          do K=J+1,nSp
            if(vCharge(K)/=0) then
              Tmp= Tmp + vMolal(J)*vMolal(K)*tZeta(J,K,jN)
            end if
          end do
        end if
      end do
      if(iDebug>2) write(ff,'(A,G15.6)') " ternary=", Tmp
      !--</ terme ternaire
      !
      AC= AC + Tmp
      !
      vLnGam(jN)=AC
    end if
  
  end do
  !--==============================</ CALCUL POUR LES ESPECES NEUTRES ==
  !
  do I=1,nSp
    if(iDebug>2) &
    & write(ff,'(I3,1X,A15,G15.6)') I,vSpc(I)%NamSp,exp(vLnGam(I))
  end do
  !
  if(iDebug>2) close(ff)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Pitzer_Calc_Gamma"
  !
  return
    
contains

  integer function IonType(I,ISIGN)
    integer,intent(in):: I      !species index
    integer,intent(in):: ISIGN  !iSign in [-1,0,1]
    
    integer:: K
    
    K= vIonType(I) != -1/0/1 for (-)/(0)/(+) respectively
    select case(ISIGN) !                !(-)  (0) (+)
      case(-1)  ; IonType= (K-1) *K/2   ! 1    0   0
      case( 0)  ; IonType= (1-K) *(1+K) ! 0    1   0
      case( 1)  ; IonType= (1+K) *K/2   ! 0    0   1
    end select
    
    return
  end function IonType

end subroutine Pitzer_Calc_Gamma

real(dp) function FONC1(X)
!--
!-- Fonction utilisee dans le calcul des coefficients d'activite
!-- par le modele de Pitzer
!--
  real(dp)::X
  if (X==Zero) then  ;  FONC1= Zero
  else               ;  FONC1= 2.0d0*(One - (One + X)*exp(-X)) /X/X
  end if
end function FONC1

real(dp) function FONC2(X)
!--
!-- Fonction utilisee dans le calcul des coefficients d'activite
!-- par le modele de Pitzer
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
!-- Calcul de la fonction ETheta(I) ==
!-- (effet electrostatique asymetrique) ==
!-- par la theorie de Pitzer (1975) ==
!-- IN ==
!  Zi: charge du premier ion
!  Zj: charge du deuxieme ion
!  SqrtI: racine carree de la force ionique
!  Aphi:  terme de Debye-Huckel
!-- OUT ==
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
!-- Element pour le calcul de la fonction ETheta(I)
!-- (effet electrostatique asymetrique)
!-- par la theorie de Pitzer (1975)
!-- IN=
!--   X : charge de l'ion considere
!-- OUT=
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
      Z=     4.0d0 * x**0.2d0    - 2.0d0
      dZ=    0.8d0 * x**(-0.8d0)
      bk(i)= Z*bk(i+1) - bk(i+2) + aik(i)
    else
      Z=     40.0d0/9.0d0 *x**(-0.1d0) - 22.0d0/9.0d0
      dZ=   -40.0d0/90.0d0*x**(-1.1d0)
      bk(i)= Z*bk(i+1) - bk(i+2) + aiik(i)
    end if
    dk(i)=  bk(i+1) + Z*dk(i+1) - dk(i+2)
  
  end do
  
  JX=  0.25d0*x - One + 0.5d0*(bk(0) - bk(2))
  J1X= 0.25d0   +    0.5d0*dZ*(dk(0) - dk(2))
  
  return
end subroutine EThetaParam

subroutine Pitzer_Calc_APhi(RhoW,EpsW,TdgK,APhi)
!--
!-- Aphi calculation is slightly different from that used in DH
!-- refs: Anderson-Crerar, p449; or Lin-Lee,FlPhasEqu-205-p69 ...
!--
  use M_Numeric_Const,   only: Pi
  use M_Dtb_Const,only: Electron,Boltzmann,Avogadro,EpsVacuum
  !
  real(dp),intent(in) :: RhoW,EpsW,TdgK
  real(dp),intent(out):: APhi
  
  Aphi= SQRT( 2*pi*Avogadro*1.D-3 ) &
  &   * (Electron/SQRT(Boltzmann) /SQRT(4.0D0*pi*EpsVacuum) )**3 &
  &   /3.0D0  &
  &   *SQRT(RhoW) &
  &   /SQRT((EpsW*TdgK)**3) &
  &   *1.0D+3
  
  !! APhi_DH= 1.824829238D6 *SQRT(RhoW) /SQRT((EpsW*TdgK)**3)
  !! call Pitzer_Calc_APhiMonnin(TdgK,1.0D0,Aphi_)
  !! print '(A,3G15.6)',"Aphi, AphiMonnin,APhi_DH:", Aphi,Aphi_,APhi_DH
  !! pause
  
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
  A0=  3.36901531d-1                 &
  &  - 6.32100430d-4 * T             &
  &  + 9.14252359    / T             & 
  &  - 1.35143986d-2 * log(T)        & 
  &  + 2.26089488d-3 / (T-263.0d0)   &
  &  + 1.92118597d-6 * T*T           &
  &  + 4.52586464d1  / (680.d0-T)     
  !
  if (T<=373.15d0) then
    !
    AV=  8.106377         &
    & - 1.256008d-1 *T    &
    & + 7.760276d-4 *T**2 &
    & - 2.098163d-6 *T**3 &
    & + 2.25777d-9  *T**4 !cm3.kg^0.5.mol(-1.5), Monnin,1989
    !
    AV=  AV * 1.0d-06 !conversion en m3
    !
  else
    !
    AV= 3.849971d2               &
    & - 6.982754     * T         &
    & + 3.877068d-2  * T*T       &
    & - 1.11381d-4   * T*T*T     &
    & + 1.589736d-7  * T*T*T*T   &
    & - 9.395266d-11 * T*T*T*T*T &
    & + 6.469241d4 / (680.d0-T) !cm3.kg^0.5.mol(-1.5), Monnin,1989
    !
    AV=  AV * 1.0d-06 !conversion en m3
    !
  end if
  !  
  Psat= 7.3649D1 -7.2582D3/T -7.3037D0*log(T) +4.1653D-6*T*T
  Psat= exp(Psat)
  !
  Aphi_ =  A0 + (Pbar*1.0D+5 - Psat) * (-AV/(4.0*8.314*T))  !Av=-1/4RT * dAphi/dP_, Monnin, 1989
  !  
  return
end subroutine Pitzer_Calc_APhiMonnin

end module M_Solmodel_Calc_Pitzer
