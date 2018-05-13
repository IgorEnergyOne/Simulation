#
import math
import numpy as np
import vectors
import space
import planet
import unit
from random import random
from constants import G
import matplotlib.pyplot as plt
import pylab

def deviation(dvx,dvy, dist):
    sp = space.Space()

    count = 0
    planetmass = 1.0 / G
    pl = planet.Planet(planetmass)
    sp.planet = pl

    unitmass = 1.0
    ut = unit.Unit(np.array([0.9,0.0,0.0]),np.array([0.0+dvx,1.0+dvy,0.0]), unitmass)
    sp.unit = ut

    t = 0

    angle = 0
    dt = 0.01 # simulation step
    pos0 = np.array([1,0.0,0.0])
    pos = ut.pos

    E0 = ut.getEnergy(planetmass)
    A0 = ut.vectorLenc(planetmass)
    U0 = ut.getPotentialEnergy(planetmass)
    KO = ut.getKinectEnergy()

#----------------Arrays of physical quantities ----------------
    arrE = []
    arrt = []
    arrT = []
    arrX = []
    arrY = []
    arrU = []
    arrK = []
    arrA = []
    arrM = []
    arrM_full = []
    arrMx = []
    arrMy = []
    arrMz = []
    arrAx = []
    arrAy = []
    arrAz = []
    arrAngle = []
    arrA_full = []
#------------------Simulation cycle-----------------------------
    n = 10 # number of circles
    while (t<(2*math.pi)*n): #one circle
        sp.step_rot(dt, dist, angle)
        count +=1
        t+=dt
        angle+= math.degrees(0.001/57)
        E = ut.getEnergy(planetmass)
        U = ut.getPotentialEnergy(planetmass)
        K = ut.getKinectEnergy()
        A = ut.vectorLenc(planetmass)
        M = ut.getMomentum()
        Angle = math.atan2(ut.vectorLenc(planetmass)[0],ut.vectorLenc(planetmass)[1])
        T = t / (2 * math.pi) # number of orbits done
#-----------------Filling the arrays----------------------------
        arrK.append(K)
        arrE.append(E)
        arrU.append(U)
        arrt.append(t)
        arrT.append(T)
        arrX.append(ut.pos[0])
        arrY.append(ut.pos[1])
        arrA.append(A)
        arrAx.append(A[0])
        arrAy.append(A[1])
        arrAz.append(A[2])
        arrA_full.append(A[0]+A[1]+A[2])
        arrMx.append(M[0])
        arrMy.append(M[1])
        arrMz.append(M[2])
        arrM_full.append(M[0]+M[1]+M[2])
        arrAngle.append(Angle)
#--------------------------------------------------------------
        #plt.plot(lnh, lnE)
#------------------	graph plotting-----------------------------
    plt.figure(1)
    plt.xlabel('Number of circles')
    plt.ylabel('Potential Energy')
    plt.title('Potential Energy of the spacecraft')
    plt.axis([0, n, 0, 2])
    plt.grid(True)
    plt.plot(arrT,arrU)
    #plt.plot(arrX,arrY)

    plt.figure(2)
    plt.xlabel('Number of circles')
    plt.ylabel('Angle')
    plt.grid(True)
    plt.plot(arrT, arrAngle)

    plt.figure(3)
    plt.xlabel('Number of circles')
    plt.ylabel('E')
    plt.title('Total Energy of the spacecraft')
    plt.axis([0, n, -0.7, -0.3])
    plt.grid(True)
    plt.plot(arrT, arrE)

    plt.figure(4)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Orbit of the spacecraft')
    plt.axis([-1.5, 1.5, -1.5, 1.5])
    plt.grid(True)
    plt.plot(arrX, arrY)

    plt.figure(5)
    plt.axis([0, n, -1, 2])
    plt.grid(True)
    plt.plot(arrT, arrA_full)
    plt.xlabel('Number of circles')
    plt.ylabel('A')
    plt.title('Laplace–Runge–Lenz vector state ')

    plt.figure(6)
    plt.axis([0, n, -1, 2])
    plt.grid(True)
    fig, ax = plt.subplots()
    ax.plot(arrT, arrAx, label ='x_component')
    ax.plot(arrT, arrAy, label ='y_component')
    ax.plot(arrT, arrAz, label ='z_component')
    plt.xlabel('Number of circles')
    plt.ylabel('A')
    plt.title('Laplace–Runge–Lenz vector components state ')
    legend = ax.legend(loc='upper right', shadow=True)

    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    plt.figure(7)
    plt.axis([0, n, -1, 2])
    plt.grid(True)
    fig, ay = plt.subplots()
    ay.plot(arrT, arrMx, label='x_component')
    ay.plot(arrT, arrMy, label='y_component')
    ay.plot(arrT, arrMz, label='z_component')
    plt.xlabel('Number of circles')
    plt.ylabel('A')
    plt.title('Angular momentum  state ')
    legend = ay.legend(loc='upper right', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('large')

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width

    plt.show()
#--------------------------------------------------------------

    M = ut.getMomentum()


    #delta = ((E - E0)**2) + vectors.squarelength(M - M0) + vectors.squarelength(A - A0)
    delta = ((E - E0) ** 2) + vectors.squarelength(A - A0)

    return delta

deviation(0,0.0,0.3)


