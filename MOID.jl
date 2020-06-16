# ======================================================================
# Description
# ========================================================================
# Calculates MOID between two keplerian orbits in AU
#
# Created: 6/16/20
# Original Author : Anivid Pedros Faura, anivid.pedrosfaura@colorado.edu
# Translated from .m to .jl by: Luke Bury, luke.bury@colorado.edu

# ========================================================================
# Libraries, Functions, and Necessary Paths
# ========================================================================
# # Libraries
# using DifferentialEquations # For using ODEProblem
# using LinearAlgebra # Allows for norm(), cross(), dot(), etc
# using Plots
using Printf # allows use of println(@sprintf())


# ========================================================================
# Dependent functions
# ========================================================================
function RefFrame(A,B)
    ### Reference frame rotation
    # 3-1-3 Tranformation (orbital elements)
    theta1 = A["Omega"]
    theta2 = A["i"]
    theta3 = A["argp"]
    RA = zeros(3,3) #from reference to A frame
    RA[1,1] = cos(theta3)*cos(theta1)-sin(theta3)*sin(theta1)*cos(theta2)
    RA[1,2] = sin(theta1)*cos(theta3)+cos(theta1)*cos(theta2)*sin(theta3)
    RA[1,3] = sin(theta2)*sin(theta3)
    RA[2,1] = -cos(theta1)*sin(theta3) - sin(theta1)*cos(theta2)*cos(theta3)
    RA[2,2] = -sin(theta3)*sin(theta1) + cos(theta3)*cos(theta2)*cos(theta1)
    RA[2,3] = sin(theta2)*cos(theta3)
    RA[3,1] = sin(theta2)*sin(theta1)
    RA[3,2] = -sin(theta2)*cos(theta1)
    RA[3,3] = cos(theta2)

    ### Preparing the orbit
    x = [cos(B["Omega"])*cos(B["argp"]) - sin(B["argp"])*cos(B["i"])*sin(B["Omega"])
         sin(B["Omega"])*cos(B["argp"]) + sin(B["argp"])*cos(B["i"])*cos(B["Omega"])
         sin(B["i"])*sin(B["argp"])]

    y = [-cos(B["Omega"])*sin(B["argp"]) - cos(B["argp"])*cos(B["i"])*sin(B["Omega"])
         -sin(B["Omega"])*sin(B["argp"]) + cos(B["argp"])*cos(B["i"])*cos(B["Omega"])
         sin(B["i"])*cos(B["argp"])];

    z = [sin(B["i"])*sin(B["Omega"])
        -sin(B["i"])*cos(B["Omega"])
        cos(B["i"])]

    xn = RA * x
    yn = RA * y
    zn = RA * z

    incliB = atan(sqrt(zn[1]^2 + zn[2]^2),zn[3])
    omegaB = -atan(zn[1],-zn[2])
    argpB = -atan(xn[3],yn[3])

    return incliB, omegaB, argpB
end


function WaterProcedure(incliB, omegaB, argpB, N, vtrueB, vL, vdis)
    if N<2
        vtrueB = zeros(4,1)
        vL = zeros(4,1)
        vdis = zeros(4,1)
        N = 4
        for j = 1:4
            vtrueB[j] = (0.25+0.5*j)*pi #evenly distributed 0.75 - 1.25 - 1.75 - 2.25
            a = (cos(omegaB)*cos(argpB+vtrueB[j]) - sin(omegaB)*sin(argpB+vtrueB[j])*cos(incliB))
            b = (sin(omegaB)*cos(argpB+vtrueB[j]) + cos(omegaB)*sin(argpB+vtrueB[j])*cos(incliB))
            vL[j] = atan(b,a)
            vdis[j] = 1e6; #set to something large
        end
    end

    return vtrueB, vL, vdis, N
end


# ========================================================================
# Calculating MOID
# ========================================================================
# -------------------------------------------------
# Setting up initial conditions
# -------------------------------------------------
# --------------------------
# Orbital parameters of body A
# --------------------------
A = Dict("sma" => 2.7691652, "e" => 0.0760091, "argp" => 73.59764*π/180, "Omega" => 80.30553*π/180, "i" => 10.59407*π/180)

# --------------------------
# Orbital parameters of body B
# --------------------------
B = Dict("sma" => 2.3655722, "e" => 0.127581, "argp" => 87.42605*π/180, "Omega" => 307.46872*π/180, "i" => 2.09575*π/180)

# -------------------------------------------------
# Calculating MOID
# -------------------------------------------------
# --------------------------
# New orbital elements: A frame
# --------------------------
# incliA = 0, omegaA = 0, argpA = 0
incliB, omegaB, argpB = RefFrame(A,B)

# --------------------------
# Scanning orbits
# --------------------------
# Scanning one full revolution of meridional plane to find local minima
# First guess for MOID is set to be big ~1e6AU

# Angle sweeping steps
cstep = 0.12 #rad - based on Wisnioski paper ~ 0.12rad
stepini = 0.07 #rad - for initial tuning step
stepfin = 1e-5 #rad - for final step of first tuning
stepmin = 1e-14 #rad - threshold step for secondtuning

# Initial Guess
trueB = -2 * cstep
moid = 1e6
dist_o = 1e6
vdis = ones(4,1)*1e6

for i in 1:2
    # First triplet
    rB = B["sma"] * (1-B["e"].^2) / (1+B["e"]*cos(trueB)) #compute the radius for B
    xB = rB * (cos(omegaB)*cos(argpB+trueB) - sin(omegaB)*sin(argpB+trueB)*cos(incliB))
    yB = rB * (sin(omegaB)*cos(argpB+trueB) + cos(omegaB)*sin(argpB+trueB)*cos(incliB))
    zB = rB * sin(argpB+trueB)*sin(incliB)

    rhoB = sqrt(xB^2+yB^2)
    L = atan(yB,xB)

    rA = A["sma"] * (1-A["e"].^2) / (1+A["e"]*cos(L)) #compute the radius for A
    rA2 = A["sma"] * (1-A["e"].^2) / (1-A["e"]*cos(L))

    if abs(rhoB-rA)>abs(rhoB+rA)
        rA = rA2
        L = L - pi
        diff = rhoB + rA2
    else
        diff = rhoB - rA
    end

    global D = zB^2+diff^2 # square of the distance
    # storing
    if i == 1
        global dist_oo = D
    else
        global dist_o = D
        global trueB_o = trueB
        global L_o = L
    end

    global trueB = trueB + cstep
end # for i in 1:2


# Full revolution
N = 0 # number of minima
dmin = D

vtrueB = fill(1,0).*0.0
vL = fill(1,0).*0.0
vdis = fill(1,0).*0.0

while trueB < (2*pi + cstep)
    rB = B["sma"] * (1-B["e"].^2) / (1+B["e"]*cos(trueB)) #compute the radius for B
    xB = rB * (cos(omegaB)*cos(argpB+trueB) - sin(omegaB)*sin(argpB+trueB)*cos(incliB))
    yB = rB * (sin(omegaB)*cos(argpB+trueB) + cos(omegaB)*sin(argpB+trueB)*cos(incliB))
    zB = rB * sin(argpB+trueB)*sin(incliB)

    rhoB = sqrt(xB^2+yB^2)
    L = atan(yB,xB)

    rA = A["sma"] * (1-A["e"].^2) / (1+A["e"]*cos(L)) #compute the radius for A
    rA2 = A["sma"] * (1-A["e"].^2) / (1-A["e"]*cos(L))

    if abs(rhoB-rA)>abs(rhoB+rA)
        rA = rA2
        L = L - pi
        diff = rhoB + rA2
    else
        diff = rhoB - rA
    end

    D = zB^2+diff^2 # square of the distance

    if (dist_o<=D) && (dist_o<=dist_oo)
        # global N = N + 1
        # vtrueB = trueB_o
        # vL = L_o
        # vdis = dist_o

        global vtrueB = append!(vtrueB,trueB_o)
        global vL = append!(vL, L_o)
        global vdis = append!(vdis, dist_o)
    end
    if dmin>D
        global dmin = D
    end

    global dist_oo = dist_o
    global trueB_o = trueB
    global L_o = L
    global dist_o = D
    global trueB = trueB + cstep
end

# -------------------------------------------------
# Water Procedure
# -------------------------------------------------
vtrueB, vL, vdis, N = WaterProcedure(incliB, omegaB, argpB, N, vtrueB, vL, vdis)

# -------------------------------------------------
# PARALLEL TUNING
# -------------------------------------------------
# Move objects separately along their orbits
# Smallest possible distance no longer meridional distance
rBt = fill(1,3).*NaN
rAt = fill(1,3).*NaN
xBt = fill(1,3).*NaN
yBt = fill(1,3).*NaN
zBt = fill(1,3).*NaN
xAt = fill(1,3).*NaN
yAt = fill(1,3).*NaN
k = 1
while k < N+2
    if k<=N
        global moid = vdis[k]
        global trueB_m = vtrueB[k]
        global L_m = vL[k]
        global myStep = stepini
        global threshold = stepfin #maybe here problem
    else
        if N == 2
            # in case of two minima are very close to each other(<1E-4 a.u.)
            # go to "water procedure"
            if (abs(vdis[1]-vdis[2])<1e-4)
                global N = 1
                global vtrueB, vL, vdis, N = WaterProcedure(incliB, omegaB, argpB, N)
                global k = 1
            else
                if (vdis[1]<moid)
                    global moid = vdis[1]
                    global trueB_m = vtrueB[1]
                    global L_m = vL[1]
                end
            end
        else
            # final tuning
            for i in 1:(N-1)
                if vdis[i]<moid
                    global moid = vdis[i]
                    global trueB_m = vtrueB[i]
                    global L_m = vL[i]
                end
            end
        end
        global myStep  = 2*stepini #inital state
        global threshold = stepmin #terminal state
    end

    rBt[2] = B["sma"] * (1-B["e"]^2) / (1+B["e"]*cos(trueB_m)) #compute the radius for B
    xBt[2] = rBt[2] * (cos(omegaB)*cos(argpB+trueB_m) - sin(omegaB)*sin(argpB+trueB_m)*cos(incliB))
    yBt[2] = rBt[2] * (sin(omegaB)*cos(argpB+trueB_m) + cos(omegaB)*sin(argpB+trueB_m)*cos(incliB))
    zBt[2] = rBt[2] * sin(argpB+trueB_m)*sin(incliB)

    rAt[2] = A["sma"] * (1-A["e"].^2) / (1+A["e"]*cos(L_m)) #compute the radius for A
    xAt[2] = rAt[2]*cos(L_m)
    yAt[2] = rAt[2]*sin(L_m)

    aleft = true
    aright = true
    bleft = true
    bright = true

    while (myStep>=threshold)
        lpoints = 0
        k1min   = 1
        k1max   = 3
        i1min   = 1
        i1max   = 3
        calc1   = false
        calc2   = false
        calc3   = false
        calc4   = false

        if bleft
            rBt[1] = B["sma"] * (1-B["e"]^2) / (1+B["e"]*cos(trueB_m-myStep)) #compute the radius for B
            xBt[1] = rBt[1] * (cos(omegaB)*cos(argpB+trueB_m-myStep) - sin(omegaB)*sin(argpB+trueB_m-myStep)*cos(incliB))
            yBt[1] = rBt[1] * (sin(omegaB)*cos(argpB+trueB_m-myStep) + cos(omegaB)*sin(argpB+trueB_m-myStep)*cos(incliB))
            zBt[1] = rBt[1] * sin(argpB+trueB_m-myStep)*sin(incliB)
            lpoints = lpoints + 1;
        end
        if bright
            rBt[3] = B["sma"] * (1-B["e"]^2) / (1+B["e"]*cos(trueB_m+myStep)) #compute the radius for B
            xBt[3] = rBt[3] * (cos(omegaB)*cos(argpB+trueB_m+myStep) - sin(omegaB)*sin(argpB+trueB_m+myStep)*cos(incliB))
            yBt[3] = rBt[3] * (sin(omegaB)*cos(argpB+trueB_m+myStep) + cos(omegaB)*sin(argpB+trueB_m+myStep)*cos(incliB))
            zBt[3] = rBt[3] * sin(argpB+trueB_m+myStep)*sin(incliB)
            lpoints = lpoints + 1;
        end
        if aleft
            rAt[1] = A["sma"] * (1-A["e"]^2) / (1+A["e"]*cos(L_m-myStep)) #compute the radius for A
            xAt[1] = rAt[1]*cos(L_m-myStep)
            yAt[1] = rAt[1]*sin(L_m-myStep)
            lpoints = lpoints + 1
        end
        if aright
            rAt[3] = A["sma"] * (1-A["e"]^2) / (1+A["e"]*cos(L_m+myStep)) #compute the radius for A
            xAt[3] = rAt[3]*cos(L_m+myStep)
            yAt[3] = rAt[3]*sin(L_m+myStep)
            lpoints = lpoints + 1
        end

        k1_t = 2
        i1_t = 2

        if lpoints == 1
            if aleft
                i1max = 1
            end
            if aright
                i1min = 3
            end
            if bright
                k1min = 3
            end
            if bleft
                k1max = 1
            end
        end

        if lpoints == 2
            if (aleft && bright)
                calc1 = true
            end
            if (aleft && bleft)
                calc2 = true
            end
            if (aright && bright)
                calc3 = true
            end
            if (aright && bleft)
                calc4 = true
            end
        end

        for k1 in k1min:k1max
            for i1 in i1min:i1max
                compute = true
                if lpoints == 2
                    if i1 != 1
                        if (k1 != 3 && calc1)||(k1 != 1 && calc2)
                            compute = false
                        end
                    end
                    if i1 != 3
                        if (k1 != 3 && calc3)||(k1 != 1 && calc4)
                            compute = false
                        end
                    end
                end
                if (k1 == 2 && i1 == 2)
                    compute = false
                end
                if (compute)
                    Dx = xBt[k1]-xAt[i1]
                    Dy = yBt[k1]-yAt[i1]
                    Dz = zBt[k1]
                    dist = Dx*Dx + Dy*Dy + Dz*Dz
                    if (dist<moid)
                        moid = dist
                        k1_t = k1
                        i1_t = i1
                    end
                end
            end
        end

        if (k1_t!=2) || (i1_t!=2)
            aleft=false
            aright=false
            bleft=false
            bright=false

            if i1_t!=2
                if i1_t == 1
                    aleft = true
                    L_m = L_m - myStep
                    rAt[3] = rAt[2]
                    xAt[3] = xAt[2]
                    yAt[3] = yAt[2]
                    rAt[2] = rAt[1]
                    xAt[2] = xAt[1]
                    yAt[2] = yAt[1]
                else
                    aright = true
                    L_m = L_m + myStep
                    rAt[1] = rAt[2]
                    xAt[1] = xAt[2]
                    yAt[1] = yAt[2]
                    rAt[2] = rAt[3]
                    xAt[2] = xAt[3]
                    yAt[2] = yAt[3]
                end
            end
            if k1_t!=2
                if k1_t == 1
                   bleft = true
                    trueB_m = trueB_m - myStep
                    rBt[3] = rBt[2]
                    xBt[3] = xBt[2]
                    yBt[3] = yBt[2]
                    zBt[3] = zBt[2]
                    rBt[2] = rBt[1]
                    xBt[2] = xBt[1]
                    yBt[2] = yBt[1]
                    zBt[2] = zBt[1]
                else
                    bright = true
                    trueB_m = trueB_m + myStep
                    rBt[1] = rBt[2]
                    xBt[1] = xBt[2]
                    yBt[1] = yBt[2]
                    zBt[1] = zBt[2]
                    rBt[2] = rBt[3]
                    xBt[2] = xBt[3]
                    yBt[2] = yBt[3]
                    zBt[2] = zBt[3]
                end
            end
        else
            aleft = true
            aright = true
            bleft = true
            bright = true
            myStep = myStep*0.15
        end
    end

    if k <= N
        vdis[k] = moid
        vtrueB[k] = trueB_m
        vL[k] = L_m
    end
    global k = k+1
end

moid = sqrt(moid)

println(@sprintf("\nThe MOID is: %1.16f\n", moid))





# ========================================================================
# Useful Formatting Structures
# ========================================================================
# -------------------------------------------------
#
# -------------------------------------------------

# --------------------------
#
# --------------------------
