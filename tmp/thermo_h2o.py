import math as m
import glpk

#help(glpk)

debug= 1

TdgK= 500. + 273.15 #temperature en Kelvin
Pbar= 1000.         #pression en bars

T0=   298.15 #temperature 'a l'e'tat T,P de reference
P0=   1.0    #pression 'a l'e'tat T,P de reference

#SUBROUTINE CalcGH2O_HolPow(TdgK,Pbar,Grt) 
#!Holland and Powell (1990), from older CdeCapitani !
#REAL(dp),INTENT(IN) :: TdgK,Pbar
#REAL(dp),INTENT(OUT):: Grt
#REAL(dp) &
#& A1,A2,A3,A4,A5,    &
#& B1,B2,B3,B4,B5,B6, &
#& C1,C2,C3,C4,C5,C6, &
#& T,T0,TT,sqT,Pkb,Pk2,Pk3,Sqp, &
#& RTLNF,             &
#& K1,K2,K3,K4,       &
#& CPRDT,CPRTDT,      &
#& H0,S0,TT0,sqT0

R_jK= 8.314510E0

H0=-241.81E0 ; S0= 188.80E-3
K1= 0.0401E0 ; K2= 0.8656E-5    ; K3=487.5E0; K4=-0.2512E0
T0= 298.15E0 ; TT0= T0*T0 #88893.4225E0

#DATA A1,A2,A3,A4,A5 &
A1= -40.338  ; A2=  1.6474    ;   A3= -0.0062115
A4=  2.0068  ; A5=  0.0562929 
#
B1=  0.117372   ; B2=  0.  ;  B3=  0.
B4= -0.00046710 ; B5=  0.  ;  B6=  0.

C1= -7.3681E-6  ;  C2= 1.10295E-7  ;  C3= -9.8774E-7
C4= -2.4819E-5  ;  C5= 8.2948E-6   ;  C6=  8.33667E-8
#
T=    TdgK
sqT0= m.sqrt(T0)
TT=   T*T
sqT=  m.sqrt(T)
#
CPRDT= K1 *(T    -T0    )      \
     + K2 *(TT   -TT0   ) /2.0 \
     - K3 *(1.0/T-1.0/T0)      \
     + K4 *(sqT  -sqT0  ) *2.0
CPRTDT= K1*m.log(T/T0)                \
      + K2*(T       -T0       )     \
      - K3*(1.0/(TT)-1.0/(TT0))/2.0 \
      - K4*(1.0/sqT -1.0/sqT0 )*2.0
#
Pkb=  Pbar/1.0E3
Pk2=  Pkb*Pkb
Pk3=  Pk2*Pkb
sqP=  m.sqrt(Pkb)
#
RTLnF=    A1 +A2*Pkb +A3*Pk2 +A4 /Pkb + A5 /Pk2                 \
        +(B1 +B2*Pkb +B3 /Pkb + B4 /Pk2  + B5 /sqP + B6 /Pk3)*T \
        +(C1 +C2*Pkb +C5 /Pkb + C3 /Pk2  + C4 /sqP + C6 /Pk3)*TT

print "RTLnF(kJ)=",RTLnF

H_atP0= H0 + CPRDT
S_atP0= S0 + CPRTDT

print "H_atP0=", H_atP0 *1.0E3
print "S_atP0=", S_atP0 *1.0E3
print "H - TS=", H_atP0 -T*S_atP0 +RTLnF*1.0E3

G= (H0 + CPRDT -T*(S0 + CPRTDT) +RTLnF) *1.0E3
G= H_atP0 -T*S_atP0 + RTLnF*1.0E3
Grt= G /R_jK/TdgK

print "G=    ",G

'''
(A) H2O Fugacities
   RT ln f(H2O) kJ/mol
   kbar\C     200    300    400    500    600    700    800    900   1000
    0.5     11.48  20.93  30.10  37.36  43.35  49.04  54.57  60.02  65.40
    1.0     12.48  22.06  31.44  39.52  46.83  53.60  60.06  66.34  72.50
    1.5     13.45  23.14  32.66  41.03  48.84  56.20  63.23  70.04  76.69
    2.0     14.40  24.18  33.81  42.38  50.46  58.17  65.57  72.75  79.75
    2.5     15.34  25.20  34.92  43.64  51.92  59.86  67.53  74.98  82.25
    3.0     16.26  26.19  35.99  44.83  53.26  61.39  69.26  76.91  84.40
    3.5     17.17  27.16  37.03  45.97  54.53  62.80  70.83  78.65  86.31
    4.0     18.06  28.11  38.04  47.07  55.73  64.12  72.29  80.26  88.06
    4.5     18.93  29.03  39.02  48.13  56.89  65.39  73.66  81.76  89.69
    5.0     19.80  29.95  39.99  49.17  58.01  66.60  74.98  83.17  91.22
'''
