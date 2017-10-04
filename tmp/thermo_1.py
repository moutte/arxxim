# compute G,H,S,V,Cp for given mineral with an EoS ala Berman
# calcul basique pour un mine'ral

import math as m

debug= 1

T= 1273.15 #temperature en Kelvin
P= 1000.   #pression en bars

T0=   298.15 #temperature 'a l'e'tat T,P de reference
P0=   1.0    #pression 'a l'e'tat T,P de reference

TT0=  T0*T0
SQT0= m.sqrt(T0)
TT=   T*T
SQT=  m.sqrt(T)

# donne'es extraites d'une base de donne'es thermodynamiques
# andalusite, Berman, 200612, 
# & ST           0.00    -2589972.10        91.4337          5.147
# & C1      236.47818      -1102.941   -7526810.000    936442368.0
# !!               K1             K4             K3             K8
# & V1     2.34426462     0.00007189    -0.07700408     0.00019235

Href= -2589972.10  # J/mole
Sref=  91.4337     # J/K/mole
Vref=  51.47       # cm3/mole
Vref=  Vref/10.    # conversion to J/bar

# coeff's for volume function
VTA=  2.34426462E-5 *Vref
VTB=  0.00007189E-5 *Vref
VPA= -0.07700408E-5 *Vref
VPB=  0.00019235E-8 *Vref

# heat capacity coeff's
K1=   236.47818  # **0
K2= 0
K3= - 7526810   # **-2
K4= - 1102.941  # **-1/2
K5= 0
K6= 0
K7= 0
K8=   936442368 # **-3
K9= 0

# heat capacity at standard T,P
FCP0 =    K1          \
        + K2 *T0      \
        + K3 /TT0     \
        + K4 /SQT0    \
        + K5 *TT0     \
        + K6 /T0      \
        + K7 *SQT0    \
        + K8 /(T0**3) \
        + K9 *T0**3

if debug: print 'FCP0=',FCP0

# heat capacity at current T,P
CPR  =    K1     \
        + K2*T   \
        + K3/TT  \
        + K4/SQT \
        + K5*TT  \
        + K6/T   \
        + K7*SQT \
        + K8/(TT*T) \
        + K9*TT*T
if debug: print 'CPR=',CPR

# integration from T0 to T (Cp dT)
CPRDT=    K1*(T     -T0     )         \
        + K2*(TT    -TT0    )/2.0     \
        - K3*(1.0/T -1.0/T0 )         \
        + K4*2*(SQT -SQT0   )         \
        + K5*(TT*T  -TT0*T0 )/3.0     \
        + K6*m.log(T/T0)              \
        + K7*(T*SQT -T0*SQT0)*2.0/3.0 \
        - K8*(1.0/TT-1.0/TT0)/2.0     \
        + K9*(TT*TT -TT0*TT0)/4.0
if debug: print 'CPRDT=',CPRDT

# integration from T0 to T (Cp/T dT)
CPRTDT=   K1*m.log(T/T0) \
        + K2*(T-T0) \
        - K3*(1.0/TT-1.0/TT0)/2.0 \
        - K4*2*(1.0/SQT-1.0/SQT0) \
        + K5*(TT-TT0)/2.0 \
        - K6*(1.0/T-1.0/T0) \
        + K7*2*(SQT-SQT0) \
        - K8*(1.0/(TT*T)-1.0/(TT0*T0))/3.0 \
        + K9*(TT*T-TT0*T0)/3.0
if debug: print 'CPRTDT=',CPRTDT

# volume computation, int from T0 to T, and from P0 to P 
FVVOL=    VTA*(T-T0)  + VTB*(T-T0)**2  \
        + VPA*(P-P0)  + VPB*(P-P0)**2
VOLUM = Vref + FVVOL

if debug: print 'VOLUM=',VOLUM

# H(T,P)= H(T,P0) + Int<P0_P>(@H@P)dP 
#       = H(T,P0) + Int<P0_P><V.dP> + T.Int<P0_P><(@V@T)dP>
# S(T,P)= S(T,P0) + Sum<P0_P>(@S@P)dP 
#       = S(T,P0)                   -   Int<P0_P><(@V@T)dP>
# G(T,P)= G(T,P0) + Int<P0_P>(@G@P)dP = G(T,P0) + Int<P0_P><V.dP>

# Int<P0_P><V.dP>
Integr_VdP= (P-P0)*(Vref +VTA*(T-T0) +VTB*(T-T0)**2)     \
          + VPA*(P**2/2.0 -P*P0    + P0**2/2.0 )          \
          + VPB*(P**3/3.0 -P**2*P0 + P*P0**2 -P0**3/3.0)

# Int<P0_P><(@V@T)dP>
Integr_dVdTdP= (VTA + 2*VTB*(T-T0)) *(P-P0)

#FCP0
#CPR

dVdT=   VTA + 2*VTB*(T - T0)
d2VdT2= 2*VTB
dVdP=   VPA + 2*VPB*(P - P0)
d2VdP2= 2*VPB

H= Href + CPRDT  # + T*Integr_dVdTdP + Integr_VdP
S= Sref + CPRTDT # - Integr_dVdTdP
G= H - T*S
G= G + Integr_VdP
H= H + Integr_VdP + T*Integr_dVdTdP 
S= S              -   Integr_dVdTdP

if debug:
  print 'results:'
  print 'H=', H
  print 'S=', S
  print 'V=', VOLUM
  print 'G=', G
  print ''
  print 'results for 1000C /1000 bars on web'
  print '(ctserver.ofm-research.org/ThermoDataSets)'
  print 'H= -2409282  J/mol'
  print 'S=  340.2731 J/K-mol'
  print 'V=  5.264202 '
  print 'G= -2842501  J/mol'

