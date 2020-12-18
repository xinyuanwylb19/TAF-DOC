# Process-based Terrestrial-aquatic DOC fluxes model (TAF-DOC)
# Version 1
# Xinyuan Wei
# 2020.10.02

import pandas as pd
import numpy as np
import math
import random
from random import randint

# Read the diver data.
wetland=pd.read_csv('Driver_P_Wetland.csv')
openwater=pd.read_csv('Driver_P_Openwater.csv')
WRT=pd.read_csv('Driver_WRT.csv')
SOC_w=pd.read_csv('Driver_SOC.csv')
WRT_R=pd.read_csv('Driver_WRT_R.csv')

temp=pd.read_csv('Driver_Temp.csv')
prep=pd.read_csv('Driver_Prep.csv')
S_dep=pd.read_csv('Driver_S_Dep.csv')
N_dep=pd.read_csv('Driver_N_Dep.csv')

# The total number of years.
T_year=34

# Outout results
DOC_Export=[]
DOC_Sediment=[]
DOC_Outgassing=[]
DOCE_Export=[]
DOCE_Sediment=[]
DOCE_Outgassing=[]

T_watersheds=len(wetland)
print ('The total number of watershed is '+str(T_watersheds)+'.')
print ('The total number of simulated years is 34.')

adj=0.1 
####################
# WSDOCM, Watershed Soil DOC Module
# The DOC is g/m2/year
####################
def WSDOCM(SOC):
    SDOC=0.1356*0.956*SOC
    return(SDOC)

####################
# WFDOCM, including five methods
# The DOC is g/m2/year
####################
def WFDOCM1(T, P, S, N):
    # Regesssion model
    if T!=0 and S!=0 and N!=0:
        y=-0.0072*T+0.0047*1.3*P-0.4039*S+0.3055*N-0.12073
        
    # Only pricipitation
    else:
        y=0.0053*P-0.1976
        
    # Adjustment for extremely low or negative value
    rd1=randint(150,200)/100
    if y<0 or y>15:
        y=0.0053*1.3*P-0.1976
        
    if y<0:
        y=0.087*rd1
        
    return (y)
    # Adjustment for extremely high value
'''
    while y>0.04:
        y=y*0.98
        if y<0.04:
            break    
'''

def WFDOCM2(T, P, S, N):
    if T!=0 and S!=0 and N!=0: 
        y=-0.0507*T+0.0256*1.3*P-0.2503*S-0.1859*N-0.6080
    else:
        y=0.0257*P-1.4293
        
    rd2=randint(150,200)/100
    if y<0 or y>15:
        y=0.0257*1.3*P-1.4293
        
    if y<0 or y>15:
        y=1.596*rd2
        
    return(y)
'''    
    while y>5:
        y=y*0.98
        if y<5:
            break
'''

def WFDOCM3(T, P, S, N):
    if T!=0 and S!=0 and N!=0:
        y=0.0906*T+0.0395*1.3*P-0.3544*S+1.7240*N-4.0303
    else:
        y=0.0369*1.3*P-2.0598
        
    rd3=randint(150,200)/100
    if y<0:
        y=0.0369*P-2.0598
        
    if y<0 or y>15:
        y=2.250*rd3
    
    return(y)
'''
    while y>7:
        y=y*0.98
        if y<7:
            break
''' 

def WFDOCM4(T, P, S, N):
    if T!=0 and S!=0 and N!=0:
        y=0.01031*T+0.1365*1.3*P-2.3890*S-0.3321*N-5.2635
    else:
        y=0.1486*1.3*P-7.7067
        
    rd4=randint(150,200)/100
    if y<0:
        y=0.1486*P-7.7067

    if y<0 or y>15:
        y=4.024*rd4
        
    return(y)
    
'''
    while y>11:
        y=y*0.98
        if y<11:
            break
'''

def WFDOCM5(T, P, S, N):
    if T!=0 and S!=0 and N!=0:
        y=-0.9553*T+0.2593*1.3*P-0.3368*S+1.3467*N-5.2082
    else:
        y=0.2238*1.3*P-20.4120
        
    rd5=randint(100,150)/100
    if y<0:
        y=0.2238*P-20.4120

    if y<0 or y>15:
        y=7.795*rd5

    return(y)
'''
    while y>13:
        y=y*0.98
        if y<13:
            break
'''

# DOC flux in the watershed. 
def WFDOC(DOC_E, p_openwater, WRT_t):
    S_rate=0.1396*math.exp(0.0121*p_openwater)
    G_rate=0.4563*WRT_t
    
    E_rate=1-S_rate-G_rate
    
    if E_rate<0.33:
        E_rate=random.randint(3000,5000)/10000
        G_rate=random.randint(2000,4000)/10000
        S_rate=1-E_rate-G_rate
        
    T_DOC=DOC_E/E_rate
    DOC_G=T_DOC*G_rate
    DOC_S=T_DOC*S_rate
         
    return(DOC_G, DOC_S)

####################
# FODOCM, DOC Fluxes to the Ocean Module
####################
def FODOCM(DOC_E,WRT_tr):
    DOCE_G=3.441*WRT_tr*DOC_E
    DOCE_S=0.0005*math.exp(148.61*WRT_tr)*DOC_E
    
    if (DOCE_G+DOCE_S)>0.32:
        DOCE_G=0.12
        DOCE_S=0.10        
    return(DOCE_G,DOCE_S)

# Estimation for each watershed.
for i in range (T_watersheds):
    print('Calculation for the watershed: '+str(i+1))
    watershed_ID=i
    p_wetland=wetland.at[watershed_ID,'P of Wetland']
    p_openwater=openwater.at[watershed_ID,'P of Openwater']
    
    WRT_t=WRT.at[watershed_ID,'WRT']
    WRT_tr=WRT_R.at[watershed_ID,'WRT_r']
    
    SOC=SOC_w.at[watershed_ID,'SOC']
    
    # Soil DOC
    SDOC=WSDOCM(SOC)

    # Decide the module
    if p_wetland<0.01:
        DOC_EM=WFDOCM1

    if p_wetland>=0.01 and p_wetland<0.05:
        DOC_EM=WFDOCM2
        
    if p_wetland>=0.05 and p_wetland<0.45:
        DOC_EM=WFDOCM3

    if p_wetland>=0.45 and p_wetland<0.55:
        DOC_EM=WFDOCM4

    if p_wetland>=0.55:
        DOC_EM=WFDOCM5

    temp_DOCG=[]
    temp_DOCS=[]
    temp_DOCE=[]
    temp_DOCEG=[]
    temp_DOCES=[]
    temp_DOCEE=[]
    
    for j in range (T_year):
        year=str(j+1985)
        Tj=float(temp.at[i,year])
        Pj=float(prep.at[i,year])
        Sj=float(S_dep.at[i,year])
        Nj=float(N_dep.at[i,year])
        
        DOC_E=DOC_EM(Tj, Pj, Sj, Nj)
        DOC_G=WFDOC(DOC_E, p_openwater, WRT_t)[0]
        DOC_S=WFDOC(DOC_E, p_openwater, WRT_t)[1]
        
        DOCE_G=FODOCM(DOC_E,WRT_tr)[0]*DOC_E
        DOCE_S=FODOCM(DOC_E,WRT_tr)[1]*DOC_E
        DOCE_E=DOC_E-DOCE_G-DOCE_S
        
        DOC_F=DOC_E+DOC_G+DOC_G
        
        # DOC_T should less than SDOC
        if DOC_F>SDOC:   
            AS_rate=0.1396*math.exp(0.0121*p_openwater)
            AG_rate=0.4563*WRT_t   
            AE_rate=1-AS_rate-AG_rate
            
            if AE_rate<0:
                AS_rate=0.21
                AG_rate=0.34
                AE_rate=0.45
            
            DOC_G=SDOC*AG_rate
            DOC_S=SDOC*AS_rate
            DOC_E=SDOC*AE_rate
            
            ADOCE_G=3.441*WRT_tr
            ADOCE_S=0.0005*math.exp(148.61*WRT_tr)
            ADOCE_E=1-ADOCE_G-ADOCE_S
            
            if ADOCE_E<0.78:
                ADOCE_G=0.12
                ADOCE_S=0.10
                ADOCE_E=0.78
        
            DOCE_G=DOC_E*ADOCE_G
            DOCE_S=DOC_E*ADOCE_S
            DOCE_E=DOC_E*ADOCE_E
        
        temp_DOCG.append(DOC_G)
        temp_DOCS.append(DOC_S)
        temp_DOCE.append(DOC_E)
        temp_DOCEG.append(DOCE_G)
        temp_DOCES.append(DOCE_S)
        temp_DOCEE.append(DOCE_E)
        
    DOC_Outgassing.append(temp_DOCG)
    DOC_Sediment.append(temp_DOCS)
    DOC_Export.append(temp_DOCE)
    DOCE_Outgassing.append(temp_DOCEG)
    DOCE_Sediment.append(temp_DOCES)
    DOCE_Export.append(temp_DOCEE)

np.savetxt('DOC_Outgassing.csv', DOC_Outgassing, delimiter=',')
np.savetxt('DOC_Sediment.csv', DOC_Sediment, delimiter=',')
np.savetxt('DOC_Export.csv', DOC_Export, delimiter=',')
np.savetxt('DOCE_Outgassing.csv', DOCE_Outgassing, delimiter=',')
np.savetxt('DOCE_Sediment.csv', DOCE_Sediment, delimiter=',')
np.savetxt('DOCE_Export.csv', DOCE_Export, delimiter=',')
