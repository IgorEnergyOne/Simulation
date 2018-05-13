#Function for simulation and algorithm of orbit optimization

import math
import numpy as np
import vectors
import space
import planet
import unit
import random
from constants import G
import matplotlib.pyplot as plt

#-----------------------------Arrays of physical quantities------------------------------
arrX = []
arrY = []
dt = 0.01  #simulation time(standart 0.001)
#----------------------------------------------------------------------------------------

def plotter(arrX, arrY, type): #function to plot graphs
    if type == 1:
        plt.figure(1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Orbit of the spacecraft')
        plt.axis([-3, 3, -3, 3])
        plt.grid(True)
        plt.plot(arrX, arrY)
        plt.show()
    if type == 2:
        plt.figure(2)
        plt.xlabel('dv_max')
        plt.ylabel('delta (10^-3)')
        plt.title('Composed function dependence on available velocity')
        plt.grid(True)
        plt.plot(arrX, arrY)
        plt.show()




def deviation(dvx1: object, dvy1: object, dvx2: object, dvy2: object, dist, plot: object = False) -> object: #simulation of first half of the orbit
    if plot:
        arrX.clear()
        arrY.clear()

    global dt
    unitmass = 1.0
    ut = unit.Unit(np.array([1.0, 0.0, 0.0]), np.array([0.0 + dvx1, 1.0 + dvy1, 0.0]), unitmass)
    sp = space.Space()

    planetmass = 1.0 / G
    pl = planet.Planet(planetmass)
    sp.planet = pl
    sp.unit = ut
    pos = ut.pos
    angle = 0

    E0 = ut.getEnergy(planetmass)
    A0 = ut.vectorLenc(planetmass)

    flag = False

    while True: #Condition for the whole orbit
        if (ut.pos[1] < ut.vel[1] * dt) and (ut.pos[0] < 0) and not flag:
            flag = True
            ut.vel[0] += dvx2
            ut.vel[1] += dvy2

        elif flag and (ut.pos[0]>0) and (ut.pos[1] > 0):
            #pos_next = ut.pos[1] + ut.vel[1] * dt
            #vel_next = ut.vel + pl.acceleration(ut.pos) * dt
            #ut.vel = (ut.vel * abs(pos_next) + vel_next * abs(ut.pos[1])) / (abs(ut.pos[1]) + abs(pos_next))
            break

        angle += math.degrees(0.001 / 57)
        sp.step_rot(dt, dist,angle)

        if plot:

            arrX.append(ut.pos[0])
            arrY.append(ut.pos[1])


    E = ut.getEnergy(planetmass)
    A = ut.vectorLenc(planetmass)
    glob_vel = ut.vel
    glob_pos = ut.pos

    delta = ((E - E0) ** 2) + vectors.squarelength(A - A0)
    return delta

#----------------------------------------------------------------------------------------------------

def deviationOpt(): # Orbit optimization
    global dt
    dist = 0.2  # distance between masses
    dvx1 = 0.0  #   initial values for fluctuations
    dvy1 = 0.0  #
    dvx2 = 0.0  #
    dvy2 = 0.0  #
    dv_opt = (0,0,0,0)
    #------------------
    delta0 = deviation(0.0,0.0,0.0,0.0, dist) #calculation of delta with no dv
    delta_min = delta0

    n = 10

    dv_range = [0, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.04,  0.05, 0.1,  0.2, 0.3]  # amount of available correction velocity
    #dv_range = [0, 0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.15]
    #step_range = [1, 10,     10,    10,    50,    50,    50,     50,    70,    80,   100,   150,   200,  300,   600,  1000, 1000, 2000]
    step_range = [1, 70, 70, 70, 80, 90, 100, 100, 150, 160, 200, 250, 300, 400, 500, 1000, 1000, 2000]
    #step_range = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
    i = -1

    print('delta0=', delta0)

    arrDelta = []
    arrVel = []
    arrVelCopm = []

    try:
        for dv in dv_range:
            count = 0
            i = i + 1
#---------------------------------------------------------------------------
            dvx1 = 0.0  # initial values for fluctuations
            dvy1 = 0.0  #
            dvx2 = 0.0  #
            dvy2 = 0.0  #
            dvx1_opt = 0.0
            dvy1_opt = 0.0
            dvx2_opt = 0.0
            dvy2_opt = 0.0
#-----------------------------------------------------------------------------
            while (count < step_range[i]):  # optimization cycle
                count += 1
                deriv = derivatives(dvx1,dvy1,dvx2,dvy2, dist)
                delta = deviation(dvx1, dvy1, dvx2, dvy2, dist, True)
                t = - delta0 / sum((i ** 2 for i in deriv))/n  # optimization step
                n = 10

#-------------------------------------------------------------------------------
#                while (opt > 10**(-3)):
#                    t = random.uniform(0, 1)
#                    opt_dvx1 = abs(deviation(dvx1-t*deriv[0], dvy1,dvx2,dvy2))
#                    opt_dvy1 = abs(deviation(dvx1, dvy1 - t * deriv[1], dvx2, dvy2))
#                    opt_dvx2 = abs(deviation(dvx1, dvy1, dvx2 - t * deriv[2], dvy2))
#                    opt_dvy2 = abs(deviation(dvx1, dvy1, dvx2, dvy2 - t * deriv[3]))
#                    opt = (opt_dvx1 + opt_dvy1 + opt_dvx2 + opt_dvy2)/4
#
#
#                    print("searching...")
#                    print("opt = "), opt
#--------------------------------------------------------------------------------

#---------------Step configuration--------------------------------------
                if (delta < 1*10**(-3)):
                    n = 50
                if (delta < 1*10**(-5)):
                    n = 1*10**2
                if (delta < 1*10**(-7)):
                    n = 1*10**4
                if (delta < 5*10**(-9)):
                    n = 1*10**6
                if (delta < 2*10**(-9)):
                    n = 1**10**7
                if (delta < 1*10**(-10)):
                    n = 5**10**7
#-----------------------------------------------------------------------

                if delta < delta_min: # in case delta in this step < delta_min
                    delta_min = delta
                    dvx1_opt = dvx1     #|
                    dvy1_opt = dvy1     #|  Save dv values to use it in simulation
                    dvx2_opt = dvx2     #|
                    dvy2_opt = dvy2     #|
                    dv_opt = (dvx1_opt, dvy1_opt, dvx2_opt, dvy2_opt)

                #arrDelta.append(math.log1p(delta))
                #arrVel.append(math.log1p(abs(dvx1) + abs(dvy1) + abs(dvx2) + abs(dvy2)))


                dvx1 += deriv[0] * t
                dvy1 += deriv[1] * t
                dvx2 += deriv[2] * t
                dvy2 += deriv[3] * t

#--------------- Projection if dv bigger than dv_max--------------------------------------
                if ((dvx1 ** 2) + (dvy1 ** 2) > (dv ** 2)):
                    K = dv / (math.sqrt((dvx1 ** 2) + (dvy1 ** 2)))
                    dvx1 = K * dvx1
                    dvy1 = K * dvy1

                if ((dvx2 ** 2) + (dvy2 ** 2) > (dv ** 2)):
                    K = dv / (math.sqrt((dvx2 ** 2) + (dvy2 ** 2)))
                    dvx2 = K * dvx2
                    dvy2 = K * dvy2
#-----------------------------------------------------------------------------------------

                print ('dvx1 = ', dvx1, ' dvy1= ', dvy1, 'dvx2=', dvx2, 'dvy2=', dvy2, ' delta = ',  delta, '   ', count)
                print ('dv1 = ', abs(dvx1)+ abs(dvy1), 'dv2 = ', abs(dvx2)+ abs(dvy2), 'dv_max= ', dv)
                print ('')

#-----------------adding result valuens to arrays-----------------------------------------
            arrDelta.append(delta_min*10**3)
            arrVel.append(dv)

            #for k in range (len(dv_range)):
            arrVelCopm.append([float(j) for j in dv_opt])
            print("arrVelComp = ", arrVelCopm)
            print('')

#-----------------------------------------------------------------------------------------
    except KeyboardInterrupt:
        print('Worked!')

    print ("delta: " ,arrDelta)
    print ("dv_max: ", arrVel)
#-----------------------------------------------------------------------------------------

    dist_range = [0.15,0.18,0.2,0.22, 0.24] #Distance for simulated systems with given dv
    arrDelta_real = []

    print("arrVelCopm length =", len(arrVelCopm))


    for dist in dist_range:
        delta_real = [0]
        arrDelta_real = [0]
        for k in range(len(arrVelCopm)-1):
            delta_real = deviation(arrVelCopm[k][0], arrVelCopm[k][1], arrVelCopm[k][2], arrVelCopm[k][3], dist)
            arrDelta_real.append(delta_real*10**3)

        plt.figure(3)
        plt.xlabel('dv_max')
        plt.ylabel('delta (10^-3)')
        plt.title('Composed function dependence on available velocity and ')
        plt.grid(True)
        plt.plot(arrVel, arrDelta_real, label = "d = 0. ")

        plt.legend()


    plt.show()


    #plotter(arrX,arrY,1)
    #plotter(arrVel, arrDelta, 2)

def derivatives(dvx1,dvy1,dvx2,dvy2, dist):
    offset = 1*(10**(-10)) #amount of offset(deviation)
    dev_nooff = deviation(dvx1,dvy1,dvx2,dvy2, dist) #derivative with no offset
    part_dx1 = (deviation(dvx1+offset,dvy1,dvx2,dvy2,dist) - dev_nooff) / offset
    part_dy1 = (deviation(dvx1, dvy1+ offset, dvx2, dvy2,dist) - dev_nooff) / offset
    part_dx2 = (deviation(dvx1,dvy1,dvx2 + offset,dvy2,dist) - dev_nooff) / offset
    part_dy2 = (deviation(dvx1, dvy1, dvx2, dvy2 + offset,dist) - dev_nooff) / offset
    return (part_dx1, part_dy1, part_dx2,  part_dy2)

