#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 20:41:18 2020

@author: chlo
"""

import numpy as np 
import spiceypy as sp

mu_sun = 1.32712440018e+11

# assuming that object B is located at the periapsis 

def DCM313(a1,a2,a3):
    c1 = np.cos(a1)
    c2 = np.cos(a2)
    c3 = np.cos(a3)
    
    s1 = np.sin(a1)
    s2 = np.sin(a2)
    s3 = np.sin(a3)
    
    d11 = c3*c1-s3*c2*s1
    d12 = c1*s1+s3*c2*c1
    d13 = s3*s2
    d21 = -s3*c1-c3*c2*s1
    d22 = -s3*s1+c3*c2*c1
    d23 = c3*s2
    d31 = s2*s1
    d32 = -s2*c1
    d33 = c2
    
    DCM = np.array([[d11, d12, d13],[d21, d22, d23],[d31, d32, d33]])
    return DCM

# water procedure - if u dont have enough minima

# used for the case in which u detect only one minimum. maybe ur steps are too big or the geometry is weird
# this is a way to force the problem to have more minimnums . if u only have 1 then this will force u to find four

def water(incB,lanB,arpB,vtrueB,vL,vdis):
    
    C0B = np.cos(lanB)
    S0B = np.sin(lanB)
    CI  = np.cos(incB)
    SI  = np.sin(incB)
    CS  = np.cos(arpB+vtrueB)
    SS  = np.sin(arpB+vtrueB)
    
    vtrueB = np.zeros([4,1]) #unless i need to store it ima not do it this way
    vL = np.zeros([4,1])
    vdis = np.zeros([4,1])
    N = np.array([0,1,2,3])
    for j in N:
        i = int(j)
        vtrueB[0][i] = (0.25+0.5*i)*np.pi #rt
        vtrueB_l= (0.25+0.5*i)*np.pi
        CA = np.cos(arpB + vtrueB)
        SA = np.sin(arpB + vtrueB)
        a = (C0B * CA) - (S0B * SA * CI)
        b = (S0B * CA) + (C0B * SA * CI)
        vL[0][i] = np.arctan2(b,a)
        vdis[0][i] = 1e6
    return np.array([vtrueB,vL,vdis,N])

    

def MOID(spk1,spk2):
    
    #get orb elems from files.. from the first state
    spiceID1r=sp.spkobj(spk1)
    for c in spiceID1r:
        spiceID1=c
        coverage=sp.spkcov(spk1,c) 
        time_interval=sp.wnfetd(coverage,0)
    cellstr1 = str(spiceID1)
    state1,lt = sp.spkezr(cellstr1,time_interval[0],'J2000','None','Sun')
    oe1=sp.oscelt(state1,time_interval[0], mu_sun)
    rp1  = oe1[0]
    ecc1 = oe1[1]
    sma1 = rp1/(1-ecc1)
    inc1 = oe1[2] #inclination
    lan1 = oe1[3] #longitude of ascending node
    arp1 = oe1[4]  # argument of periapsis
    
    spiceID2r=sp.spkobj(spk2)
    for c in spiceID2r:
        spiceID2=c
        coverage2=sp.spkcov(spk2,c) 
        time_interval2=sp.wnfetd(coverage2,0)
    cellstr2 = str(spiceID2)
    state2,lt = sp.spkezr(cellstr2,time_interval2[0],'J2000','None','Sun')
    oe2=sp.oscelt(state2,time_interval2[0], mu_sun)
    rp2  = oe2[0]
    ecc2 = oe2[1]
    sma2 = rp2/(1-ecc2)
    inc2_i = oe2[2] #inclination
    lan2_i = oe2[3] #longitude of ascending node
    arp2_i = oe2[4]  # argument of periapsis

    # redefining A orbit frame as reference frame by representing A frames wrt inertial, and then recalculating B orbitOEs
    
    #first step - defining A and B frames wrt inertial using a 3-1-3 EA rotation
    dcm1 = DCM313(inc1,lan1,arp1) 
    dcm2_i = DCM313(inc2_i,lan2_i,arp2_i) 
    
    #inertial is just [0 0 0 ; 0 0 0 ; 0 0 0], thus, these DCMs are also the base frame vectors for the A,B frame wrt inertial,
    #note that dcm1 is the DCM for translating vectors from inertial frame to ref frame A
    
    #second step - recalculating the correspnding OEs
    
    
    #this is the base frame vectors for B frame wrt A frame (reference)
    dcm2_x = dcm1 @ np.array([[dcm2_i[0][0]],[dcm2_i[1][0]],[dcm2_i[2][0]]])
    dcm2_y = dcm1 @ np.array([[dcm2_i[0][1]],[dcm2_i[1][1]],[dcm2_i[2][1]]])
    dcm2_z = dcm1 @ np.array([[dcm2_i[0][2]],[dcm2_i[1][2]],[dcm2_i[2][2]]])
    
    dcm2 = np.array([dcm2_x,dcm2_y,dcm2_z])
    
    #getting the new inc, lan, arp for B ; Schaub&Jenkins p94 eq 3.37
    # this is diff from anivids - check why
    
    inc2 = np.arctan2(np.sqrt(((dcm2[2][0])**2 + (dcm2[2][1])**2)),(dcm2[2][2])) #arccos c33
    lan2 = -np.arctan2(dcm2[2][0],-dcm2[2][1])
    arp2 = -np.arctan2((dcm2[0][2]),(dcm2[1][2]))
    
    # angle sweeping steps 
    cstep = 0.12
    
    # generating a good initial guess .. initializing the dist making sure there r no problems
    trueB = -2 * cstep #true anomaly of B, negative as if it were past info s.t. can start running true anom at 0 later
    dist_o = 1e6 #initial distance btwn points before time step
    vdis = np.ones([4,1]) * 1e6 
    
    listi = np.array([1,2])
    
    C0B = np.cos(lan2)
    S0B = np.sin(lan2)
    CI  = np.cos(inc2)
    SI  = np.sin(inc2)
    CS  = np.cos(arp2+trueB)
    SS  = np.sin(arp2+trueB)
    
    for i in listi:
        print ('first bit',i)
        
#        C0B = np.cos(lan2)
#        S0B = np.sin(lan2)
#        CI  = np.cos(inc2)
#        SI  = np.sin(inc2)
#        CS  = np.cos(arp2+trueB)
#        SS  = np.sin(arp2+trueB)
        
        rB = sma2 * (1-ecc2**2)/(1+(ecc2*np.cos(trueB))) #computing radius
        xB = rB * (C0B*CS-S0B*SS*CI) #from paper
        yB = rB * (S0B*CS+C0B*SS*CI)
        zB = rB * (SS*SI)
        
        rhoB = np.sqrt(xB**2 + yB**2)
        L = np.arctan(yB, xB) #this is the longitude of A, this is how A is rotating as the true anom of B changes
        
        rA = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L))  #this is rA if L is ok
        rA2 = sma1 * (1-ecc1**2)/(1-ecc1*np.cos(L)) #this is rA if L is not ok. cos then L-pi works
        
        # this is checking for the case where obj A is at the opposite side of the Sun
        # basically verifying that your guess of obj A location is ok
        # if its not, then you just do L-pi to move body A to the opposite/mirror pt of the ellipse
        if abs(rhoB-rA)>abs(rhoB+rA):
            rA   = rA2
            L    = L-np.pi
            diff = rhoB + rA2
        else:
            diff = rhoB - rA
        
        D = zB ** 2 + diff**2
        
        #this is where u generate guesses to initalize the next part (the while loop)
        # d_o will be the distance btwn the two point on orbit A and B, one time step before
        # d_oo will be the distance two time steps before
        if i == 1:
            dist_oo = D
        else:
            dist_o = D #this is just defining the first distance ur working w .. at the very first iteration i=0..
            trueB_o = trueB #this is just defining the first true anom u are working w
            L_o = L #defining the first longitude of A
        trueB = trueB + cstep
    
    N=0 #the number of minima found
    dmin=D #starting with this first distance
    vtrueB = np.zeros((4,1))
    vL = np.zeros((4,1))
    vdis = np.zeros((4,1))
    
    # scanning of the orbit using the min distance to start with 
    while trueB < (2*np.pi + cstep): #while we have not made the full revolution yet
        print('trueB each step in revolution while loop', trueB)
        #compute B radius
        rB = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB))
        xB = rB * (C0B*CS-S0B*SS*CI)
        yB = rB * (S0B*CS+C0B*SS*CI)
        zB = rB * (SS*SI)
        
        rhoB = np.sqrt(xB**2 + yB**2)
        L = np.arctan2(yB, xB)
        
        rA = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L))
        rA2 = sma1 * (1-ecc1**2)/(1-ecc1*np.cos(L))
        
        # this is checking for the case where obj A is at the opposite side of the Sun
        # basically verifying that your guess of obj A location is ok
        # if its not, then you just do L-pi to move body A to the opposite/mirror pt of the ellipse
        if abs(rhoB-rA)>abs(rhoB+rA):
            rA   = rA2
            L    = L-np.pi
            diff = rhoB + rA2
        else:
            diff = rhoB - rA
        
        D = zB ** 2 + diff**2
        
        # if old distance dist_o is smaller than both old old distance d_oo and recent distance D
        # then it is a local minimum. store that info.
        if dist_o <= D and dist_o <= dist_oo:
            N = N + 1
            vtrueB[N] = trueB_o
            vL[N] = L_o
            vdis[N] = dist_o
        #otherwise if dmin is bigger than D, then we set the new dmin to be D
        #lets us keep establishing what is the lowest D found so far
        if dmin>D:
            dmin=D
        # readjusting the things we are iterting, so old d becomes old old, and recent becomes old,etc
        dist_oo=dist_o
        trueB_o = trueB
        L_o=L
        dist_o = D
        #new iteration!
        trueB = trueB + cstep
        
        
        #so this first part finds the first four candiates for moid
        #getting an estim for the MOiD
        # but this is not definitely true so in the next part , for every min we found, we are trying to find what point is actually the minimum in the proximity of that point
        #NEXT PART. BY SEVERAL iterations trying to get to point a and b minimum
        # taking the point found from first part as a starting pont 
        
        #initializing water procedure if minima is too few
    if trueB == 2*np.pi: # @ last iteration
        if N < 2: #idk if this is a correct limit
            water(inc2,lan2,arp2,trueB,vL,vdis) #call upon water procedure
    
    # the coordinates of the six poinnts for each each orbit at that 
    #recall we work w 6 points, two points and then one on the left n one on the right for both
                    
    stepini = 0.07
    stepfin = 1e-5
    stepmin = 1e-14
    moid = 1e6
    
    
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
    
    k=0
    
    # moving through all the minimua
    while k < N+2: #plus two bc you hvae two additional points .. gotta loop through three points (the plus and minus right and left points)
        # trying to move slightly from the ideal plane that we first used and trying to move a bit farther and closer .. moving off the meridional plane 
        print ('k value for line 269 while loop',k)
        if k <= N: # for each minimum, we are setting the moid to be that distance..  assuming each minimum is the right answer and then testing.
            moid = vdis[k]
            print ('this is moid line 286', moid)
            trueB_m = vtrueB[k] 
            L_m = vL[k]
            step = stepini #changing the step size through
            threshold = stepfin
        else:
            if N == 2 : #just a special case
                if (np.abs(vdis[0] - vdis[1]) < 1e-4): #if the two minimums are super close
                    N = 1 #should this be zero bc python starts at 0
                    [vtrueB, vL, vdis, N] = water(inc2, lan2, arp2, N) #do water procedure again, starting over the procedure bc what u had is not right
                    k = 1 #should this be zero bc python starts at 0
                else:
                    if (vdis[0] < moid): #updating the moid if the new point found here is less than the current MOID, then update the minimum distance stored
                        moid = vdis[0]
                        trueB_m = vtrueB[0]
                        L_m = vL[0]
            else:
                #so if N is not equal to two , keeps updating as well 
                #running thru this mostly
                # storing all the MOIDS awe have and choosign minium one
                i_list = np.array[0:N-1] #does this make a list from 1 to N-1 
                for i in i_list:
                    if vdis[0] < moid:
                        moid = vdis[0]
                        trueB_m = vtrueB[i]
                        L_m = vL[i]
            # based on what is in the paper  to cover up more space 
            # this is to cover all the area
            step = 2*stepini
#            print ('step line 314',step)
            threshold = stepmin
        
        C0B = np.cos(lan2)
        S0B = np.sin(lan2)
        CI  = np.cos(inc2)
        SI  = np.sin(inc2)
        CSM = np.cos(arp2+trueB_m)
        SSM = np.sin(arp2+trueB_m)           
        
        #using the second index, [1] because we are using the point in the middle as reference
        # [0] and [2] are the right and left htings 
        rBt[1] = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB_m))
        xBt[1] = rBt[1] *  ((C0B*CSM)-(S0B*SSM*CI))
        yBt[1] = rBt[1] *  ((S0B*CSM)+(C0B*SSM*CI))
        zBt[1] = rBt[1] *  (SSM*SI)
        
        rAt[1] = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L_m))
        xAt[1] = rAt[1]*np.cos(L_m)
        yAt[1] = rAt[1]*np.sin(L_m)
        
        #gtting more precise here
        
        #saying were going to use all the points
        aleft = True
        aright = True
        bleft = True
        bright = True
        
        #while youre still relatively far form the point
        #until you end up very close to the actual point
        
        if step >= threshold: #switched while to if statement see if thisdoes anyting
#            print('step threshold comparison',step)
            lpoints = 0 #number of extra points im computing
            # refer to the B orbit , the left and right. for these ponts. if min is 1 and max is three. 
            #maybe after keep iterationg. the right direction doesnt minimize. so if aright is false so for the rest of the looping not looking at the right part
            #checking if in that direction, the distance increases. and determines whether or not going in a direction is worth it 
            #if the saddle behavior happns where u go increase an then down 
            #if the steps are small enough then technically u wouldve found a local min by finding that
            #can have weird behavior.. but it shold be rare that it happens that u miss stuff
            k1min = 1
            k1max = 3
            #refer to A orbit
            i1min = 1
            i1max = 3
            #so at first we're calculating right, left, then center
        
            #these indicate .. so u need to check left or right point yes or no
            #and that sets the min or max that i need to use
            # checking after first iteration, if aleft bleft etc.. 
            #at firs tu say they are true 

            # in the next interation after the while it will keep getting rid of points im not going to use
            calc1 = False
            calc2 = False
            calc3 = False
            calc4 = False
            
            #in these statemetn we're just obtaining the coordinates for that point
            #and with that later on i'm looking at which of these distances is the minimum
            # so we have six points six coords w each situation
            # and testing which distance between these,
            # all the combos possible. to see which combo gives the minumum
            
            #thats why subtracting or adding step. going to right or left
            
            #these are placing the points in the space
            #jut trying to put the points in the map
            if bleft:
#                print ('bleft')
                rBt[0] = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB_m-step))
                xBt[1] = rBt[1] *  ((C0B*np.cos(arp2 + trueB_m - step))-(S0B*np.sin(arp2 + trueB_m - step)*CI))
                yBt[1] = rBt[1] *  ((S0B*np.cos(arp2 + trueB_m - step))+(C0B*np.sin(arp2 + trueB_m - step)*CI))
                zBt[1] = rBt[1] *  (np.sin(arp2 + trueB_m - step)*SI)  
                lpoints = lpoints + 1 
            if bright:
#                print ('bright')
                rBt[2] = sma2 * (1-ecc2**2)/(1+ecc2*np.cos(trueB_m+step))
                xBt[2] = rBt[2] *  ((C0B*np.cos(arp2 + trueB_m + step))-(S0B*np.sin(arp2 + trueB_m + step)*CI))
                yBt[2] = rBt[2] *  ((S0B*np.cos(arp2 + trueB_m + step))+(C0B*np.sin(arp2 + trueB_m + step)*CI))
                zBt[2] = rBt[2] *  (np.sin(arp2 + trueB_m + step)*SI)
                lpoints = lpoints + 1 
            if aleft:
#                print ('aleft')
                rAt[0] = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L_m-step))
                xAt[0] = rAt[0] * np.cos(L_m-step)
                yAt[0] = rAt[0] * np.sin(L_m-step)
                lpoints = lpoints + 1
            if aright:
#                print ('aright')
                rAt[2] = sma1 * (1-ecc1**2)/(1+ecc1*np.cos(L_m+step))
                xAt[2] = rAt[2] * np.cos(L_m+step)
                yAt[2] = rAt[2] * np.sin(L_m+step)
                lpoints = lpoints + 1
        #talking about the main middle point
        #using the middle positoin as a reference
        # the index for the last candidate thats been found
        k1_t = 2
#        print ('k1_t=2')
        i1_t = 2
        
        #if in the previous iteraiton you got rid of a lot of c
        #ases
        #so if u hav
        #if i got rid of some points bc u realized that 
        #just to make sure that after the new iteration you're 
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
        forlist = np.arange(k1min,k1max) #yeaa i dont think this works like that
        
        # trying to fin which combo of points gives the minumum 
        for k1 in forlist:
#            print ('k1 in forlist line 444',k1)
            forlist2 = np.arange(i1min,i1max) #lol again. this doesnt work but
            for i1 in forlist2:
                compute = True
                #this next part checking to see if based on prev info u have u dont have to bother computing the distances 
                #instead of having in total six points u hve four points 
                # so if lpoint == 2 , not considering one of the sides for each of the ref points
                if lpoints == 2: #this means u r only considering one point on each side of each reference. so ur either just considering right to one and left to one, left ot one right to one etc
                    if i1 != 1: #this means if its NOT equal to 
                        if (k1 != 3 and calc1) or (k1 != 1 and calc2):  #whaaaaaaaat
                            compute = False
                    if i1 != 3: #this means if its NOT 
                        if (k1 != 3 and calc1) or (k1!=1 and calc4): #whaaaaaat is happening in anis code
                            compute = False
                #if its the one in the middle for A and the one in the middle fo B
                #dont compute bc its already ur reference
                if (k1 == 2 and i1 == 2):
                    compute = False
                # if u do have to compute
                #then compute it
                # k is the letter to place the diff b points. i is for A points. just doing normal distance
                # A is not going to have any z component bc of the way we framed the ref plane
                if compute:
                    Dx = xBt[k1] - xAt[i1]
                    Dy = yBt[k1] - yBt[i1]
                    Dz = zBt[k1]
                    dist = Dx*Dx+Dy*Dy+Dz*Dz
                    if dist<moid: #every timeu compute, u check to see if this new combo is lower to the reference 
                        moid = dist
                        print ('this is the moid line 473',moid)
                        k1_t = k1 # the combo that makes a minimum has this position 
                        i1_t = i1
        if k1_t !=2 or i1_t!=2: #if ur answer was not the same as the current reference point, the ones in hte middle, then update it 
#            print ('reached line 476')
            aleft = False
            aright = False
            bleft = False
            bright = False
            #for each of them, if the ref was not the right one, update so that this new thing is the new considered minumum 
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
        else: #before changing the new step, evaluate all the possibilities again , want the whole pic again starting from a diff point
            print ('reached line 502')
            aleft = True
            aright = True
            bleft = True
            bright = True
            # as long as step is bigger than min distance set threshold .. then we keep looking until we find the final point we want for that minimum
            step = step * 0.15 #for every minimums ur testing, every iteraiton, you are changing the step size .. reducing the points between the two orbtis, bc ur assuming u are converging 
        if k <= N:
            vdis[k] = moid #vdis is a vector used to define the list of MOIDS until u find the right MOID
            vtrueB[k] = trueB_m
            vL[k] = L_m
        k = k+1
    moid = np.sqrt(moid)
    print (moid)
    return np.array([moid]) #idk
    

JANUSMOIDFilesToLoad = ['200117_Janus_FG3_Open_LowThrust.bsp', 
               '200117_Janus_VH_Open_LowThrust.bsp'] 

sp.furnsh('de430.bsp')
sp.furnsh('naif0012tls.txt')
sp.furnsh(JANUSMOIDFilesToLoad) 




MOID('200117_Janus_FG3_Open_LowThrust.bsp','200117_Janus_VH_Open_LowThrust.bsp')