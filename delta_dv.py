#Function for simulation and algorithm of orbit optimization

import math
import numpy as np
import vectors
import space
import planet
import unit
from random import random
from constants import G
import matplotlib.pyplot as plt

#-----------------------------Arrays of physical quantities------------------------------
arrX = []
arrY = []
glob_pos = []
glob_vel = []
#----------------------------------------------------------------------------------------

def plotter(arrX, arrY, type): #function to plot graphs
    if type == 1:
        plt.figure(1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Orbit')
        plt.axis([-2, 2, -2, 2])
        plt.grid(True)
        plt.plot(arrX, arrY)
        plt.show()
    if type == 2:
        plt.figure(2)
        plt.xlabel('dv')
        plt.ylabel('delta')
        plt.title('Chaos graph')
        plt.grid(True)
        plt.plot(arrX, arrY)
        plt.show()


def deviation(dvx,dvy,dist): #simulation of first half of the orbit
    global glob_pos     #some global names to store variables
    global glob_vel

    unitmass = 1.0
    ut = unit.Unit(np.array([1.0, 0.0, 0.0]), np.array([0.0 + dvx, 1.1 + dvy, 0.0]), unitmass)
    sp = space.Space()

    count = 0
    planetmass = 1.0 / G
    pl = planet.Planet(planetmass)
    sp.planet = pl
    sp.unit = ut
    pos = ut.pos

    dt = 0.001 #simulation step

    E0 = ut.getEnergy(planetmass)
    A0 = ut.vectorLenc(planetmass)
    M0 = ut.getMomentum()

    while ((ut.pos[1] > ut.vel[1] * dt) or (ut.pos[0] > 0 )): #Condition for one half of the orbit
    #while (ut.movetime(dt) < math.pi):
    #while (count < math.pi/dt):
        sp.step(dt, dist)
        count +=1
        arrX.append(ut.pos[0])
        arrY.append(ut.pos[1])


    E = ut.getEnergy(planetmass)
    A = ut.vectorLenc(planetmass)
    M = ut.getMomentum()
    glob_vel = ut.vel
    glob_pos = ut.pos

    #delta = ((E - E0)**2) + vectors.squarelength(M - M0) + vectors.squarelength(A - A0)
    delta = ((E - E0) ** 2) + vectors.squarelength(A - A0)
    #print ('unit position: ',ut.pos)
    return delta


#----------------------------------------------------------------------------------------------------

def deviation2(dvx,dvy,dist): # simulation of the second half
    sp = space.Space()

    count = 0
    planetmass = 1.0 / G
    pl = planet.Planet(planetmass)
    sp.planet = pl

    unitmass = 1.0

    pos = glob_pos
    vel = glob_vel
    dv = np.array([dvx,dvy,0])

    ut = unit.Unit(pos,(vel+dv), unitmass)
    sp.unit = ut
    #print ('start velocity', ut.vel)
    dt = 0.001

    E0 = ut.getEnergy(planetmass)

    A0 = ut.vectorLenc(planetmass)
    M0 = ut.getMomentum()

    #while (ut.pos[1] < ut.vel[1] * dt) or (ut.pos[0] < 0 ): #Condition for second half of the orbit
    #while (ut.pos[1] <= 0):
    while (count < math.pi/dt):
        sp.step(dt, dist)
        count +=1
        arrX.append(ut.pos[0])
        arrY.append(ut.pos[1])

    E = ut.getEnergy(planetmass)
    A = ut.vectorLenc(planetmass)
    M = ut.getMomentum()


    #delta = ((E - E0)**2) + vectors.squarelength(M - M0) + vectors.squarelength(A - A0)
    delta = ((E - E0) ** 2) + vectors.squarelength(A - A0)
    #print ('unit position: ', ut.pos , 'unit velocity: ', ut.vel)
    return delta

#----------------------------------------------------------------------------------------------------

def deviationOpt(): # trajectroty optimization on one half of the orbit
    d = 0.01  #step of the optimization
    dist = 0.0
    dvx0 = 0.0  #initial fluctuations
    dvy0 = 0.0  #
    delta0 = deviation(0.0,0.0, dist)
    delta = 0.0
    count = 0.0
    arrDelta = []
    arrVel = []
    arrStep = np.arange(1.75*(10**(-3)),1*(10**(-3)),-1*(10**(-4)))

    print ( 'initial delta =', delta0)
    try:
        while (delta0 > 4.51 * (10 ** (-7))):  # cycle of first half optimization
            dvx = dvx0 + d * (random() * 2 - 1)
            dvy = dvy0 + d * (random() * 2 - 1)
            count += 1

            delta = deviation(dvx, dvy, dist)

            if (delta < 4.6*10**(-7)):
                d = 0.005

            if (delta < delta0):
                dvx0 = dvx
                dvy0 = dvy
                delta0 = delta

                # if (count > 75):
                # break
            print('dvx0 = ', dvx0, ' dvy0= ', dvy0, " delta = ", delta, "   ", count)

        print('optimization one complete')
        print('dvx0 = ', dvx0, ' dvy0= ', dvy0, " delta = ", delta, "   ")
        print('')

        #for deltaSet in arrStep:
        while (delta0 > 6 * 10**(-12)):  # cycle of second half optimization
            dvx = dvx0 + d * (random() * 2 - 1)
            dvy = dvy0 + d * (random() * 2 - 1)
            count += 1

            delta = deviation2(dvx, dvy, dist)

            if (delta < 3.6 * 10 ** (-10)):
                d = 0.0005

            if (delta < delta0):
                dvx0 = dvx
                dvy0 = dvy
                delta0 = delta
            print('dvx0 = ', dvx0, ' dvy0= ', dvy0, " delta = ", delta, "   ", count)

        print('optimization two complete')
        print('')
        print("dvx0 = ", dvx0, " dvy0= ", dvy0, " delta = ", delta, "   ", " delta0 = ", delta0, "   ")
    except KeyboardInterrupt:
        print('Worked!')

    print(len(arrX), len(arrY))

    plotter(arrVel,arrDelta,2)

