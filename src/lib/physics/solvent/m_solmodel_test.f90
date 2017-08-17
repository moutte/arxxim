module M_Solmodel_Test
!--
  use M_Kinds
  
  implicit none
  
  private
  
  public:: Pitzer_EThetaCalc_Test
  
contains
  
subroutine Pitzer_EThetaCalc_Test
!--
!-- test de la routine de calcul EThetaCalc,
!-- pour différentes valeurs de Zi (- 2,3,4)
!-- -> résulatst tabulés dans "etheta_test.tab"
!--
!-- Calcul de la fonction ETheta(I) --
!-- (effet electrostatique asymetrique) --
!-- par la theorie de Pitzer (1975) --
!
!-- IN --
!  Zi: charge du premier ion
!  Zj: charge du deuxieme ion
!  SqrtI: racine carree de la force ionique
!  Aphi:  terme de Debye-Huckel
!
!-- OUT --
!  Eth  : effet electrostatique asymetrique
!  Eth1 : derivee de Eth par rapport a la temperature
!--
  use M_IoTools
  use M_Solmodel_Calc_Pitzer
  !
  integer :: F,I
  real(dp):: IonicStr, APHI
  real(dp):: eth, eth1
  !
  character, parameter:: tt= char(9)
  
  call GetUnit(f)
  open(f,file="etheta_test.tab")
  
  write(f,'(8(A,A1))') &
  & ".I",char(9),"IonicStr",char(9), &
  & "-Eth",char(9),"IonicStr*Eth1",char(9), &
  & "-Eth",char(9),"IonicStr*Eth1",char(9), &
  & "-Eth",char(9),"IonicStr*Eth1",char(9)
  !
  Aphi= 0.30D0
  IonicStr= 1.D-6
  
  I= 1
  do
    
    write(f,'(I3,A1,G15.6,A1)',advance="NO") &
    & I,tt,IonicStr,tt
    
    call EThetaCalc (2, 1, SQRT(IonicStr), Aphi, Eth, Eth1)
    write(f,'(2(G15.6,A1))',advance="NO") &
    & -Eth,tt,IonicStr*Eth1,tt
    
    call EThetaCalc (3, 1, SQRT(IonicStr), Aphi, Eth, Eth1)
    write(f,'(2(G15.6,A1))',advance="NO") &
    & -Eth,tt,IonicStr*Eth1,tt
    
    call EThetaCalc (4, 1, SQRT(IonicStr), Aphi, Eth, Eth1)
    write(f,'(2(G15.6,A1))') &
    & -Eth,tt,IonicStr*Eth1,tt
    
    IonicStr= IonicStr *1.5D0
    I= I+1
    if(IonicStr>20.0D0) exit
    
  end do
  
  close(f)
  
  return
end subroutine Pitzer_EThetaCalc_Test

end module M_Solmodel_Test
