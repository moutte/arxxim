MODULE M_Solmodel_Test
!--
  USE M_Kinds
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC:: Pitzer_EThetaCalc_Test
  
CONTAINS
  
SUBROUTINE Pitzer_EThetaCalc_Test
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
  USE M_IoTools
  USE M_Solmodel_Calc_Pitzer
  !
  INTEGER :: F,I
  REAL(dp):: IonicStr, APHI
  REAL(dp):: eth, eth1
  !
  CHARACTER, PARAMETER:: tt= CHAR(9)
  
  CALL GetUnit(f)
  OPEN(f,FILE="etheta_test.tab")
  
  WRITE(f,'(8(A,A1))') &
  & ".I",CHAR(9),"IonicStr",CHAR(9), &
  & "-Eth",CHAR(9),"IonicStr*Eth1",CHAR(9), &
  & "-Eth",CHAR(9),"IonicStr*Eth1",CHAR(9), &
  & "-Eth",CHAR(9),"IonicStr*Eth1",CHAR(9)
  !
  Aphi= 0.30D0
  IonicStr= 1.D-6
  
  I= 1
  DO
    
    WRITE(f,'(I3,A1,G15.6,A1)',ADVANCE="NO") &
    & I,tt,IonicStr,tt
    
    CALL EThetaCalc (2, 1, SQRT(IonicStr), Aphi, Eth, Eth1)
    WRITE(f,'(2(G15.6,A1))',ADVANCE="NO") &
    & -Eth,tt,IonicStr*Eth1,tt
    
    CALL EThetaCalc (3, 1, SQRT(IonicStr), Aphi, Eth, Eth1)
    WRITE(f,'(2(G15.6,A1))',ADVANCE="NO") &
    & -Eth,tt,IonicStr*Eth1,tt
    
    CALL EThetaCalc (4, 1, SQRT(IonicStr), Aphi, Eth, Eth1)
    WRITE(f,'(2(G15.6,A1))') &
    & -Eth,tt,IonicStr*Eth1,tt
    
    IonicStr= IonicStr *1.5D0
    I= I+1
    IF(IonicStr>20.0D0) EXIT
    
  ENDDO
  
  CLOSE(f)
  
  RETURN
END SUBROUTINE Pitzer_EThetaCalc_Test

ENDMODULE M_Solmodel_Test
