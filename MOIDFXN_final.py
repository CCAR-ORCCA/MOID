#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 14:48:07 2020

@author: chlo
"""

import numpy as np 
import spiceypy as sp

mu_sun = 1.32712440018e+11


def DCM313(a1,a2,a3):
    c1 = np.cos(a1)
    c2 = np.cos(a2)
    c3 = np.cos(a3)
    
    s1 = np.sin(a1)
    s2 = np.sin(a2)
    s3 = np.sin(a3)
    
    d11 = c3*c1-s3*c2*s1
    d12 = c3*s1+s3*c2*c1
    d13 = s3*s2
    d21 = -s3*c1-c3*c2*s1
    d22 = -s3*s1+c3*c2*c1
    d23 = c3*s2
    d31 = s2*s1
    d32 = -s2*c1
    d33 = c2
    
    DCM = np.array([[d11, d12, d13],[d21, d22, d23],[d31, d32, d33]])
    return DCM

# water procedure 
def water(incB,lanB,arpB,vtrueB,vL,vdis):
    
    C0B = np.cos(lanB)
    S0B = np.sin(lanB)
    CI  = np.cos(incB)
    SI  = np.sin(incB)
    CS  = np.cos(arpB+vtrueB)
    SS  = np.sin(arpB+vtrueB)
    
    vtrueB = [] 
    vL = []
    vdis = []
    N = 4
    NLIST = np.array([0,1,2,3])
    for j in NLIST:
        i = int(j)
        vtrueB.append((0.25+0.5*i)*np.pi) 
        vtrueB_l= (0.25+0.5*i)*np.pi
        CA = np.cos(arpB + vtrueB)
        SA = np.sin(arpB + vtrueB)
        a = (C0B * CA) - (S0B * SA * CI)
        b = (S0B * CA) + (C0B * SA * CI)
        vL.append(np.arctan2(b,a))
        vdis.append(1e6)
    return np.array([vtrueB,vL,vdis,N])

def MOID(spk1,spk2):
#    spiceID1r=sp.spkobj(spk1)
#    for c in spiceID1r:
#        spiceID1=c
#        coverage=sp.spkcov(spk1,c) 
#        time_interval=sp.wnfetd(coverage,0)
#    cellstr1 = str(spiceID1)
#    state1,lt = sp.spkezr(cellstr1,time_interval[0],'J2000','None','Sun')
    oe1=np.array([2.7691652, 0.0760091,  10.59407,  80.30553,73.59764])
    sma1=oe1[0]
    ecc1 = oe1[1]
    inc1 = np.radians(oe1[2]) #inclination
    lan1 = np.radians(oe1[3]) #longitude of ascending node
    arp1 = np.radians(oe1[4])  # argument of periapsis
#    spiceID2r=sp.spkobj(spk2)
#    for c in spiceID2r:
#        spiceID2=c
#        coverage2=sp.spkcov(spk2,c) 
#        time_interval2=sp.wnfetd(coverage2,0)
#    cellstr2 = str(spiceID2)
#    state2,lt = sp.spkezr(cellstr2,time_interval2[0],'J2000','None','Sun')
    oe2 = np.array([2.3655722, 0.127581 ,2.09575,  307.46872, 87.42605, ])
    sma2=oe2[0]
    ecc2 = oe2[1]
    inc2_i = np.radians(oe2[2]) #inclination
    lan2_i = np.radians(oe2[3]) #longitude of ascending node
    arp2_i = np.radians(oe2[4])  # argument of periapsis

    dcm1 = DCM313(lan1,inc1,arp1) 
    dcm2_i = DCM313(lan2_i,inc2_i,arp2_i) 

    dcm2_x = dcm1 @ np.array([[dcm2_i[0][0]],[dcm2_i[0][1]],[dcm2_i[0][2]]])
    dcm2_y = dcm1 @ np.array([[dcm2_i[1][0]],[dcm2_i[1][1]],[dcm2_i[1][2]]])
    dcm2_z = dcm1 @ np.array([[dcm2_i[2][0]],[dcm2_i[2][1]],[dcm2_i[2][2]]])
    
    dcm2 = np.array([dcm2_x,dcm2_y,dcm2_z])
    
    inc2 = np.arctan2(np.sqrt(((dcm2[2][0])**2 + (dcm2[2][1])**2)),(dcm2[2][2])) 
    lan2 = -np.arctan2(dcm2[2][0],-dcm2[2][1])
    arp2 = -np.arctan2((dcm2[0][2]),(dcm2[1][2]))
    
    cstep = 0.12
    stepini = 0.07
    stepfin = 1e-5
    stepmin = 1e-14
    
    trueB = -2 * cstep 
    moid=1e6
    dist_o = 1e6 
    vdis = [1e6,1e6,1e6,1e6]
    
    listi = np.array([1,2])
    
    C0B = np.cos(lan2)
    S0B = np.sin(lan2)
    CI  = np.cos(inc2)
    SI  = np.sin(inc2)
    CS  = np.cos(arp2+trueB)
    SS  = np.sin(arp2+trueB)
    
    for i in listi:
        C0Bi = np.cos(lan2)
        S0Bi = np.sin(lan2)
        CIi  = np.cos(inc2)
        SIi  = np.sin(inc2)
        CSi  = np.cos(arp2+trueB)
        SSi  = np.sin(arp2+trueB)
        rB = sma2 * (1-ecc2**2)/(1+(ecc2*np.cos(trueB))) 
        
        xB = rB * (C0Bi*CSi-S0Bi*SSi*CIi)
        yB = rB * (S0Bi*CSi+C0Bi*SSi*CIi)
        zB = rB * (SSi*SIi)
        
        rhoB = np.sqrt(xB**2 + yB**2)
        L = np.arctan2(yB, xB)

        rA = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L))  
        rA2 = sma1 * (1-ecc1**2)/(1-ecc1*np.cos(L)) 
        
        if abs(rhoB-rA)>abs(rhoB+rA):
            rA   = rA2
            L    = L-np.pi

            diff = rhoB + rA2
        else:
            diff = rhoB - rA
        D = zB ** 2 + diff**2

        if i == 1:
            dist_oo = D
        else:
            dist_o = D 
            trueB_o = trueB 
            L_o = L
        trueB = trueB + cstep
    
    N=0 
    dmin=D
    vtrueB=[1,1,1,1]
    vL=[1,1,1,1]
    

    while trueB < ((2*np.pi) + cstep): 
        C0Bq = np.cos(lan2)
        S0Bq = np.sin(lan2)
        CIq  = np.cos(inc2)
        SIq  = np.sin(inc2)
        CSq  = np.cos(arp2+trueB)
        SSq  = np.sin(arp2+trueB)
        rB = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB))
        xB = rB * (C0Bq*CSq-S0Bq*SSq*CIq)
        yB = rB * (S0Bq*CSq+C0Bq*SSq*CIq)
        zB = rB * (SSq*SIq)
        
        rhoB = np.sqrt(xB**2 + yB**2)
        L = np.arctan2(yB, xB)
        
        rA = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L))
        rA2 = sma1 * (1-ecc1**2)/(1-ecc1*np.cos(L))

        if abs(rhoB-rA)>abs(rhoB+rA):
            rA   = rA2
            L    = L-np.pi
            diff = rhoB + rA2
           
        else:
            diff = rhoB - rA
        
        D = zB ** 2 + diff**2
        
        if dist_o <= D and dist_o <= dist_oo:
            N = N + 1
            vtrueB[N-1] = trueB_o
            vL[N-1]=L_o
            vdis[N-1]=dist_o

        if dmin>D:
            dmin=D
       
        dist_oo=dist_o
        trueB_o = trueB
        L_o=L
        dist_o = D
        
        trueB = trueB + cstep
        
    if trueB == 2*np.pi: 
        if N < 2: 
            water(inc2,lan2,arp2,trueB,vL,vdis) 
#            
    
    rBt = np.empty((3,1,))
    rAt = np.empty((3,1,))
    xBt = np.empty((3,1,))
    xAt = np.empty((3,1,))
    yBt = np.empty((3,1,))
    yAt = np.empty((3,1,))
    zBt = np.empty((3,1,))
    
    rBt[:] = np.nan
    rAt[:] = np.nan
    xBt[:] = np.nan
    xAt[:] = np.nan
    yBt[:] = np.nan
    yAt[:] = np.nan
    zBt[:] = np.nan
    
    k=1
    
    
    while k < N+2:         
        if k <= N:            
            moid = vdis[k-1]
            trueB_m = vtrueB[k-1] 
            L_m = vL[k-1]
            step = stepini 
            threshold = stepfin
        else:
            if N == 2 :
                if (np.abs(vdis[0] - vdis[1]) < 1e-4): 
                    N = 1 
                    [vtrueB, vL, vdis, N] = water(inc2, lan2, arp2, N,vL,vdis) 
                    k = 1 
                else:
                    if (vdis[0] < moid): 
                        moid = vdis[0]
                        trueB_m = vtrueB[0]
                        L_m = vL[0]
            else:
                i_list = np.arange(0,N-1) 
                for i in i_list:
                    if vdis[0] < moid:
                        moid = vdis[0]
                        trueB_m = vtrueB[i]
                        L_m = vL[i]
        
            step = 2*stepini
            threshold = stepmin
        
        C0B = np.cos(lan2)
        S0B = np.sin(lan2)
        CI  = np.cos(inc2)
        SI  = np.sin(inc2)
        CSM = np.cos(arp2+trueB_m)
        SSM = np.sin(arp2+trueB_m)           
        
        rBt[1] = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB_m))
        xBt[1] = rBt[1] *  ((C0B*CSM)-(S0B*SSM*CI))
        yBt[1] = rBt[1] *  ((S0B*CSM)+(C0B*SSM*CI))
        zBt[1] = rBt[1] *  (SSM*SI)
        rAt[1] = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L_m))
        xAt[1] = rAt[1]*np.cos(L_m)
        yAt[1] = rAt[1]*np.sin(L_m)
        
        aleft = True
        aright = True
        bleft = True
        bright = True
        
        while step >= threshold: 
            lpoints = 0 
            k1min = 1
            k1max = 3
            #refer to A orbit
            i1min = 1
            i1max = 3
            
            # in the next interation after the while it will keep getting rid of points im not going to use
            calc1 = False
            calc2 = False
            calc3 = False
            calc4 = False
            if bleft:
                rBt[0] = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB_m-step))
                xBt[0] = rBt[0] *  ((C0B*np.cos(arp2 + trueB_m - step))-(S0B*np.sin(arp2 + trueB_m - step)*CI))
                yBt[0] = rBt[0] *  ((S0B*np.cos(arp2 + trueB_m - step))+(C0B*np.sin(arp2 + trueB_m - step)*CI))
                zBt[0] = rBt[0] *  (np.sin(arp2 + trueB_m - step)*SI)  
                lpoints = lpoints + 1 
            if bright:
                rBt[2] = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB_m+step))
                xBt[2] = rBt[2] *  ((C0B*np.cos(arp2 + trueB_m + step))-(S0B*np.sin(arp2 + trueB_m + step)*CI))
                yBt[2] = rBt[2] *  ((S0B*np.cos(arp2 + trueB_m + step))+(C0B*np.sin(arp2 + trueB_m + step)*CI))
                zBt[2] = rBt[2] *  (np.sin(arp2 + trueB_m + step)*SI)
                lpoints = lpoints + 1 
            if aleft:
                rAt[0] = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L_m-step))
                xAt[0] = rAt[0] * np.cos(L_m-step)
                yAt[0] = rAt[0] * np.sin(L_m-step)
                lpoints = lpoints + 1
            if aright:
                rAt[2] = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L_m+step))
                xAt[2] = rAt[2] * np.cos(L_m+step)
                yAt[2] = rAt[2] * np.sin(L_m+step)
                lpoints = lpoints + 1
    
            k1_t = 2
            i1_t = 2
            
            if lpoints == 1:
                if aleft:
                    i1max = 1
                if aright:
                    i1min = 3
                if bright:
                    k1min = 3
                if bleft:
                    k1max = 1
                    
            if lpoints == 2:
                if aleft and bright:
                    calc1 = True
                if aleft and bleft:
                    calc2 = True
                if aright and bright:
                    calc3 = True
                if aright and bleft :
                    calc4 = True
                    
            forlist = np.arange(k1min,k1max+1) 
            forlist2 = np.arange(i1min,i1max+1)
            
            for k1 in forlist:
                for i1 in forlist2:
                    compute = True
                    if lpoints == 2: 
                        if i1 != 1: 
                            if (k1 != 3 and calc1) or (k1 != 1 and calc2): 
                                compute = False
                        if i1 != 3:
                            if (k1 != 3 and calc1) or (k1!=1 and calc4): 
                                compute = False
                    if (k1 == 2 and i1 == 2):
                        compute = False
                    if compute:
                        Dx = xBt[k1-1] - xAt[i1-1]                        
                        Dy = yBt[k1-1] - yAt[i1-1]
                        Dz = zBt[k1-1]
                        dist = Dx*Dx+Dy*Dy+Dz*Dz
                        if dist<moid: 
                            moid = dist
                            k1_t = k1 
                            i1_t = i1
            if k1_t !=2 or i1_t!=2:
                aleft = False
                aright = False
                bleft = False
                bright = False
                if i1_t != 2: 
                    if i1_t == 1:
                        aleft = True
                        L_m = L_m - step
                        rAt[2] = rAt[1]
                        xAt[2] = xAt[1]
                        yAt[2] = yAt[1]
                        rAt[1] = rAt[0]
                        xAt[1] = xAt[0]
                        yAt[1] = yAt[0]
                    else:
                        aright = True
                        L_m = L_m + step
                        rAt[0] = rAt[1]
                        xAt[0] = xAt[1]
                        yAt[0] = yAt[1]
                        rAt[1] = rAt[2]
                        xAt[1] = xAt[2]
                        yAt[1] = yAt[2]
                if k1_t != 2:
                    if k1_t == 1:
                        bleft = True
                        trueB_m = trueB_m - step
                        rBt[2] = rBt[1]
                        xBt[2] = xBt[1]
                        yBt[2] = yBt[1]
                        zBt[2] = zBt[1]
                        rBt[1] = rBt[0]
                        xBt[1] = xBt[0]
                        yBt[1] = yBt[0]
                        zBt[1] = zBt[0]
                    else:
                        bright = True
                        trueB_m = trueB_m + step
                        rBt[0] = rBt[1]
                        xBt[0] = xBt[1]
                        yBt[0] = yBt[1]
                        zBt[0] = zBt[1]
                        rBt[1] = rBt[2]
                        xBt[1] = xBt[2]
                        yBt[1] = yBt[2]
                        zBt[1] = zBt[2]
                        
                        
            else:
                aleft = True
                aright = True
                bleft = True
                bright = True
                step = step * 0.15 
            
        if k <= N:
            vdis[k-1] = moid          
            vtrueB[k-2] = trueB_m
            vL[k-2] = L_m
        k = k+1
    moid = np.sqrt(moid)
    print ('THIS IS MOID',moid)
    return [moid]
    

JANUSMOIDFilesToLoad = ['200117_Janus_FG3_Open_LowThrust.bsp', 
               '200117_Janus_VH_Open_LowThrust.bsp'] 

sp.furnsh('de430.bsp')
sp.furnsh('naif0012tls.txt')
sp.furnsh(JANUSMOIDFilesToLoad) 



MOID('200117_Janus_FG3_Open_LowThrust.bsp','200117_Janus_VH_Open_LowThrust.bsp')