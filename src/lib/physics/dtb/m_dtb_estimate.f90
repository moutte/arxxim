MODULE M_Dtb_Estimate
!--
!-- estimating thermodynamic properties of silicate minerals,
!-- especially clay minerals,
!-- by interpolation methods, e.g. Chermak - Rimstidt
!--
  USE M_Kinds
  !USE M_IoTools
  USE M_Trace,ONLY: iDebug,T_,fTrc,Stop_,Pause_
  !USE M_T_DtbAquHkf,ONLY: T_DtbAquHkf
  !USE M_T_DtbMinHkf,ONLY: T_DtbMinHkf
  !USE M_T_DtbMinThr,ONLY: T_DtbMinThr
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dtb_Estim_Phyllit
  PUBLIC:: Dtb_Estim_Init
  !  
  REAL(dp):: &
  G_X_K2O,G_X_NA2O,G_X_CAO, &
  G_O_FEO,G_O_FEOH2,G_O_MGO,G_O_MGOH2,G_O_CAO, &
  G_O_AL2O3,G_O_ALOH3,G_T_AL2O3,G_T_SIO2,G_H2O,G_O_FE2O3, &
  H_X_K2O,H_X_NA2O,H_X_CAO, &
  H_O_FEO,H_O_FEOH2,H_O_MGO,H_O_MGOH2,H_O_CAO, &
  H_O_AL2O3,H_O_ALOH3,H_T_AL2O3,H_T_SIO2,H_H2O,H_O_FE2O3
  !
CONTAINS

SUBROUTINE Chermak_Init
  !
  G_X_K2O=   -722.94D0  ;  H_X_K2O=   -735.24D0 !X_K2O     
  G_X_NA2O=  -672.50D0  ;  H_X_NA2O=  -683.00D0 !X_NA2O    
  G_X_CAO=   -710.08D0  ;  H_X_CAO=   -736.04D0 !X_CAO     
  G_O_FEO=   -266.29D0  ;  H_O_FEO=   -290.55D0 !O_FEO     
  G_O_FEOH2= -542.04D0  ;  H_O_FEOH2= -596.07D0 !O_FE(OH)2 
  G_O_MGO=   -628.86D0  ;  H_O_MGO=   -660.06D0 !O_MGO     
  G_O_MGOH2= -851.86D0  ;  H_O_MGOH2= -941.62D0 !O_MG(OH)2 
  G_O_CAO=   -669.13D0  ;  H_O_CAO=   -696.65D0 !O_CAO     
  G_O_AL2O3=-1594.52D0  ;  H_O_AL2O3=-1690.18D0 !O_AL2O3   
  G_O_ALOH3=-1181.62D0  ;  H_O_ALOH3=-1319.55D0 !O_AL(OH)3 
  G_T_AL2O3=-1631.32D0  ;  H_T_AL2O3=-1716.24D0 !T_AL2O3   
  G_T_SIO2=  -853.95D0  ;  H_T_SIO2=  -910.97D0 !T_SIO2    
  G_H2O=     -239.91D0  ;  H_H2O=     -292.37D0 !H2O       
  G_O_FE2O3= -776.07D0  ;  H_O_FE2O3= -939.18D0 !O_FE2O3   
  !
ENDSUBROUTINE

SUBROUTINE Dtb_Estim_Phyllit(vEle)
!--
!-- generate dGf,dHf for phyllites according to Chermak_Rimstidt
!--
  USE M_IoTools,  ONLY: GetUnit
  USE M_Files,    ONLY: DirOut
  USE M_T_Element,ONLY: T_Element,Formula_Build,Element_Index
  !
  TYPE(T_Element),DIMENSION(:),INTENT(IN):: vEle
  !
  !REAL(dp),DIMENSION(14):: vG,vH
  INTEGER :: iK_,iMg,iFe,iAl,iSi
  INTEGER :: I,J,Div,ZSp,F
  CHARACTER(LEN=255):: S
  INTEGER,DIMENSION(:),ALLOCATABLE:: vStoik
  !
  REAL(dp):: G0,H0
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dtb_Calc_Phyllit"
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"phyllite.res")
  !
  ALLOCATE(vStoik(1:SIZE(vEle))); vStoik=0
  !
  iK_= Element_Index("K__",vEle) ; IF(iK_==0) RETURN
  iMg= Element_Index("MG_",vEle) ; IF(iMg==0) RETURN
  iAl= Element_Index("AL_",vEle) ; IF(iAl==0) RETURN
  iSi= Element_Index("SI_",vEle) ; IF(iSi==0) RETURN
  !
  !ILLITE = K(X+Y) MG(Y)AL(2-Y)  AL(X)SI(4-X)  O(12)H(2)
  !       = K(I+J) MG(I)AL(12-I) AL(J)SI(24-J) O(72)H(12) DIV(6)
  !
  !X=0,Y=0 -> _.__AL2.SI4   -> PYROPHYLLITE
  !X=1,Y=0 -> K.__AL2.ALSI3 -> MUSCOVITE
  !X=0,Y=1 -> K.MGAL_.SI4   -> CELADONITE
  !
  Div= 6
  ZSp= 0
  DO I=0,6
    DO J=0,6
      IF(I+J<7) THEN
        !
        vStoik(iK_)= I+J     
        vStoik(iMg)= I       
        vStoik(iAl)= J-I+12  
        vStoik(iSi)= 24-J
        !
        CALL Formula_Build(vEle,vStoik,Zsp,Div,S)
        !
        !X_K2O=    REAL(I+J )/12.D0   !X_K2O=    (I+J)  /6       /2
        !O_MGO=    REAL(I   )/ 9.D0   !O_MGO=     I     /6 *(2/3)
        !O_MG(OH)= REAL(I   )/18.D0   !O_MG(OH)=  I     /6 *(1/3)
        !O_AL2O3=  REAL(12-I)/18.D0   !O_AL2O3=  (12-I) /6 *(2/3)/2
        !O_AL(OH)3=REAL(12-I)/18.D0   !O_AL(OH)3=(12-I) /6 *(1/3)
        !T_AL2O3=  REAL(J   )/12.D0   !T_AL2O3=   J     /6       /2
        !T_SIO2=   REAL(24-J)/ 6.D0   !T_SIO2=   (24 -J)/6       /2
        G0= G_X_K2O   *REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        & + G_O_MGO   *REAL(I   )/ 9.D0 &  !O_MGO=     I     /6 *(2/3)
        & + G_O_MGOH2 *REAL(I   )/18.D0 &  !O_MG(OH)=  I     /6 *(1/3)
        & + G_O_AL2O3 *REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        & + G_O_ALOH3 *REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        & + G_T_AL2O3 *REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        & + G_T_SIO2  *REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
        H0= H_X_K2O   *REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        & + H_O_MGO   *REAL(I   )/ 9.D0 &  !O_MGO=     I     /6 *(2/3)
        & + H_O_MGOH2 *REAL(I   )/18.D0 &  !O_MG(OH)=  I     /6 *(1/3)
        & + H_O_AL2O3 *REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        & + H_O_ALOH3 *REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        & + H_T_AL2O3 *REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        & + H_T_SIO2  *REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
        !!PRINT '(A,G15.6)',TRIM(S),G0  ; PAUSE_
        !!PRINT '(A,G15.6)',TRIM(S),H0  ; PAUSE_
        !
        WRITE(f,'(A,A1,2(G15.6,A1))') TRIM(S),T_, -G0/LOG(10.D0),T_, H0,T_
        !
      ENDIF
    ENDDO
  ENDDO
  !
  vStoik= 0
  iFe= Element_Index("FE_",vEle); IF(iFe==0) RETURN
  DO I=0,6
    DO J=0,6
      IF(I+J<7) THEN
        vStoik(iK_)= I+J     
        vStoik(iFe)= I       
        vStoik(iAl)= J-I+12  
        vStoik(iSi)= 24-J
        !
        CALL Formula_Build(vEle,vStoik,Zsp,Div,S)
        G0= G_X_K2O   *REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        & + G_O_FEO   *REAL(I   )/ 9.D0 &  !O_FEO=     I     /6 *(2/3)
        & + G_O_FEOH2 *REAL(I   )/18.D0 &  !O_FE(OH)=  I     /6 *(1/3)
        & + G_O_AL2O3 *REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        & + G_O_ALOH3 *REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        & + G_T_AL2O3 *REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        & + G_T_SIO2  *REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
        !!PRINT '(A,G15.6)',TRIM(S),G0  ; PAUSE_
        !
        H0= H_X_K2O   *REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        & + H_O_FEO   *REAL(I   )/ 9.D0 &  !O_MGO=     I     /6 *(2/3)
        & + H_O_FEOH2 *REAL(I   )/18.D0 &  !O_MG(OH)=  I     /6 *(1/3)
        & + H_O_AL2O3 *REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        & + H_O_ALOH3 *REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        & + H_T_AL2O3 *REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        & + H_T_SIO2  *REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
        !!PRINT '(A,G15.6)',TRIM(S),H0  ; PAUSE_
        !
        WRITE(f,'(A,A1,2(G15.6,A1))') TRIM(S),T_, -G0/LOG(10.D0),T_, H0,T_
        !
      ENDIF
    ENDDO
  ENDDO
  !
  DEALLOCATE(vStoik)
  !  
  CLOSE(F)
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dtb_Calc_Phyllit"
ENDSUBROUTINE Dtb_Estim_Phyllit

SUBROUTINE Dtb_Estim_Init(vSpc)
  USE M_T_Species,ONLY: T_Species,Species_Index
  !
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  !
  INTEGER:: I
  !INTEGER:: I_Kaoli,I_Musco,I_Talc_,I_Pyrof,I_Celad,I_Kfels,I_Chrys,I_Phlog
  INTEGER, DIMENSION(1:8):: vI
  REAL(dp),DIMENSION(1:8):: vG
  !
  INTEGER,DIMENSION(1:8,1:8):: tStoik= RESHAPE((/  &
  !  !X_K2O, X_NA2O, O_MGO  O_MG(OH) O_AL2O3  O_AL(OH)3 T_AL2O3 T_SIO2
  &  0,      0,      0,     0,       -12,     12,       0,      0,  &    !KAOL   
  &  4,      -20,    6,     0,       18,      0,        14,     -3, &    !MUSC   
  &  0,      24,     0,     0,       0,       0,        0,      0,  &    !PARAGON
  &  -8,     -8,     8,     -4,      0,       0,        8,      0,  &    !TALC   
  &  -16,    -16,    -6,    0,       6,       -6,       -2,     3,  &    !PYROPH 
  &  -4,     -4,     -6,    0,       -18,     0,        10,     3,  &    !MICROCL
  &  0,      0,      -4,    8,       0,       0,        0,      0,  &    !CHRYSOT
  &  24,     24,     0,     0,       0,       0,        -24,    0/),&    !CELADMG
  (/8,8/) ) 
  !
  vI(1)= Species_Index("KAOLINITE",   vSpc) !I_Kaoli
  vI(2)= Species_Index("MUSCOVITE",   vSpc) !I_Musco
  vI(3)= Species_Index("PARAGONITE",   vSpc) !I_
  vI(4)= Species_Index("TALC",        vSpc) !I_Talc_
  vI(5)= Species_Index("PYROPHYLLITE",vSpc) !I_Pyrof
  vI(6)= Species_Index("MICROCLINE",  vSpc) !I_Kfels
  vI(7)= Species_Index("CHRYSOTILE",  vSpc) !I_Chrys
  vI(8)= Species_Index("CELADONITE",  vSpc) !I_Celad
  !
  PRINT '(/,A,/)', "Dtb_Estim_Init"
  PRINT '(I3,2A)',vI(1),"=","KAOLINITE"
  PRINT '(I3,2A)',vI(2),"=","MUSCOVITE"
  PRINT '(I3,2A)',vI(3),"=","PARAGONITE"
  PRINT '(I3,2A)',vI(4),"=","TALC"
  PRINT '(I3,2A)',vI(5),"=","PYROPHYLLITE"
  PRINT '(I3,2A)',vI(6),"=","MICROCLINE"
  PRINT '(I3,2A)',vI(7),"=","CHRYSOTILE"
  PRINT '(I3,2A)',vI(8),"=","CELADONITE"
  CALL Pause_
  !
  
  !  !X_K2O, X_NA2O, O_MGO  O_MG(OH) O_AL2O3  O_AL(OH)3 T_AL2O3 T_SIO2
  !  0,      0,      0,     0,       -6,      6,        0,      0      !KAOL   
  !  2,      -10,    3,     0,       9,       0,        7,      -1.5   !MUSC   
  !  0,      12,     0,     0,       0,       0,        0,      0      !PARAGON
  !  -4,     -4,     4,     -2,      0,       0,        4,      0      !TALC   
  !  -8,     -8,     -3,    0,       3,       -3,       -1,     1.5    !PYROPH 
  !  -2,     -2,     -3,    0,       -9,      0,        5,      1.5    !MICROCL
  !  0,      0,      -2,    4,       0,       0,        0,      0      !CHRYSOT
  !  12,     12,     0,     0,       0,       0,        -12,    0      !CELADMG
  !the original stoichio table (multiplied by 6 !!)
  !(table is the inverse of the transpose of this matrix)
  !   __	X_K2O	X_NA2O	O_MGO	O_MG(OH	O_AL2O3 O_AL(OH	T_AL2O3 T_SIO2
  !   KAOL	0	0	0	0	2	8	0	12
  !   MUSC	3	0	0	0	4	4	3	18
  !   PARAG	0	3	0	0	4	4	3	18
  !   TALC	0	0	12	6	0	0	0	24
  !   PYROPH	0	0	0	0	4	4	0	24
  !   MICROCL	3	0	0	0	0	0	3	18
  !   CHRYSOT	0	0	6	12	0	0	0	12
  !   CELADMG	3	0	4	2	2	2	0	24

  !
  DO i= 1,8
    PRINT '(I3)',tStoik(1,i) ! /12.0D0
  ENDDO
  CALL Pause_
  DO i= 1,8
    vG(i)= DOT_PRODUCT(vSpc(vI(1:8))%G0rt,tStoik(i,1:8)) /12.0D0
  ENDDO
  G_X_K2O=   vG(1) !;  H_X_K2O=  vH(1)    !X_K2O     
  G_X_NA2O=  vG(2) !;  H_X_NA2O= vH(2)    !X_NA2O    
  G_O_MGO=   vG(3) !;  H_O_MGO=  vH(3)    !O_MGO     
  G_O_MGOH2= vG(4) !;  H_O_MGOH2=vH(4)    !O_MG(OH)2 
  G_O_AL2O3= vG(5) !;  H_O_AL2O3=vH(5)    !O_AL2O3   
  G_O_ALOH3= vG(6) !;  H_O_ALOH3=vH(6)    !O_AL(OH)3 
  G_T_AL2O3= vG(7) !;  H_T_AL2O3=vH(7)    !T_AL2O3   
  G_T_SIO2=  vG(8) !;  H_T_SIO2= vH(8)    !T_SIO2    
  !
  DO i= 1,8
    PRINT '(G15.6)',vG(i)
  ENDDO
  CALL Pause_
  !
  !G_X_CAO=         !;  H_X_CAO=      !X_CAO     
  !G_O_FEO=         !;  H_O_FEO=      !O_FEO     
  !G_O_FEOH2=       !;  H_O_FEOH2=    !O_FE(OH)2 
  !G_O_CAO=         !;  H_O_CAO=      !O_CAO     
  !G_H2O=           !;  H_H2O=        !H2O       
  !G_O_FE2O3=       !;  H_O_FE2O3=    !O_FE2O3   
  !
ENDSUBROUTINE Dtb_Estim_Init

ENDMODULE M_Dtb_Estimate

  !old! vG(1)=  -722.94D0;  vH(1)=  -735.24D0 !X_K2O 
  !old! vG(2)=  -672.50D0;  vH(2)=  -683.00D0 !X_NA2O 
  !old! vG(3)=  -710.08D0;  vH(3)=  -736.04D0 !X_CAO  
  !old! vG(4)=  -266.29D0;  vH(4)=  -290.55D0 !O_FEO    
  !old! vG(5)=  -542.04D0;  vH(5)=  -596.07D0 !O_FE(OH)2
  !old! vG(6)=  -628.86D0;  vH(6)=  -660.06D0 !O_MGO    
  !old! vG(7)=  -851.86D0;  vH(7)=  -941.62D0 !O_MG(OH)2
  !old! vG(8)=  -669.13D0;  vH(8)=  -696.65D0 !O_CAO    
  !old! vG(9)= -1594.52D0;  vH(9)= -1690.18D0 !O_AL2O3  
  !old! vG(10)=-1181.62D0;  vH(10)=-1319.55D0 !O_AL(OH)3
  !old! vG(11)=-1631.32D0;  vH(11)=-1716.24D0 !T_AL2O3  
  !old! vG(12)= -853.95D0;  vH(12)= -910.97D0 !T_SIO2   
  !old! vG(13)= -239.91D0;  vH(13)= -292.37D0 !H2O      
  !old! vG(14)= -776.07D0;  vH(14)= -939.18D0 !O_FE2O3  
        !
        !G0= vG(1 )*REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        !& + vG(6 )*REAL(I   )/ 9.D0 &  !O_MGO=     I     /6 *(2/3)
        !& + vG(7 )*REAL(I   )/18.D0 &  !O_MG(OH)=  I     /6 *(1/3)
        !& + vG(9 )*REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        !& + vG(10)*REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        !& + vG(11)*REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        !& + vG(12)*REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
        !H0= vH(1 )*REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        !& + vH(6 )*REAL(I   )/ 9.D0 &  !O_MGO=     I     /6 *(2/3)
        !& + vH(7 )*REAL(I   )/18.D0 &  !O_MG(OH)=  I     /6 *(1/3)
        !& + vH(9 )*REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        !& + vH(10)*REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        !& + vH(11)*REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        !& + vH(12)*REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
        !
        !old! G0= vG(1 )*REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        !old! & + vG(4 )*REAL(I   )/ 9.D0 &  !O_FEO=     I     /6 *(2/3)
        !old! & + vG(5 )*REAL(I   )/18.D0 &  !O_FE(OH)=  I     /6 *(1/3)
        !old! & + vG(9 )*REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        !old! & + vG(10)*REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        !old! & + vG(11)*REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        !old! & + vG(12)*REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !old! !
        !old! !!PRINT '(A,G15.6)',TRIM(S),G0  ; PAUSE_
        !old! !
        !old! H0= vH(1 )*REAL(I+J )/12.D0 &  !X_K2O=    (I+J)  /6       /2
        !old! & + vH(6 )*REAL(I   )/ 9.D0 &  !O_MGO=     I     /6 *(2/3)
        !old! & + vH(7 )*REAL(I   )/18.D0 &  !O_MG(OH)=  I     /6 *(1/3)
        !old! & + vH(9 )*REAL(12-I)/18.D0 &  !O_AL2O3=  (12-I) /6 *(2/3)/2
        !old! & + vH(10)*REAL(12-I)/18.D0 &  !O_AL(OH)3=(12-I) /6 *(1/3)
        !old! & + vH(11)*REAL(J   )/12.D0 &  !T_AL2O3=   J     /6       /2
        !old! & + vH(12)*REAL(24-J)/ 6.D0    !T_SIO2=   (24 -J)/6       /2
        !
