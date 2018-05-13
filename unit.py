# unit module
import math
import numpy as np
import vectors as v

G = 6.67408*(10.0**(-11))
EULER = 0
EULER_BETTER = 1
RUNGE_KUTTA = 2


class Unit:
    def __init__(self, pos=np.array([0.0,0.0,0.0]), vel=np.array([0.0,0.0,0.0]), mass=1.0, method = EULER):
        self.mass = mass
        self.pos = pos
        self.vel = vel
        self.method = method


    def getEnergy(self, M):
        return self.mass*((v.length(self.vel))**2)/2 - G*M*self.mass/v.length(self.pos)

    def getPotentialEnergy(self, M):
        return G*M*self.mass/v.length(self.pos)

    def getKinectEnergy(self):
        return  self.mass*((v.length(self.vel))**2)/2

    def vectorLenc(self, M):
        return np.cross(np.cross(self.vel,self.pos),self.vel)*self.mass*self.mass - self.pos*self.mass*G*M/v.length(self.pos)

    def getMomentum(self):
        return np.cross(self.pos,self.vel)*self.mass


    def move(self, a, dt, dist):
        if self.method == EULER:
            self.vel = self.vel + a(self.pos, dist) * dt
            self.pos = self.pos + self.vel * dt
        elif self.method == EULER_BETTER:
            pos_prev = self.pos
            pos_next = self.pos + self.vel * dt
            a_prev = a(pos_prev, dist)
            a_next = a(pos_next, dist)
            vel_prev = self.vel
            vel_next = self.vel + a_next * dt
            self.pos = self.pos + (vel_next + vel_prev) * dt / 2
            self.vel = self.vel + (a_next + a_prev) * dt / 2
        elif self.method == RUNGE_KUTTA:
            k1 = self.vel
            k2 = self.vel + a(self.pos+k1*(dt/2),dist)*(dt/2)
            k3 = self.vel + a(self.pos+k2*(dt/2),dist)*(dt/2)
            k4 = self.vel + a(self.pos+k3*dt, dist)*dt
            self.pos = self.pos + dt/6*(k1 + 2*k2 + 2*k3 + k4)
            self.vel = k4

    def move_rot(self, a, dt, dist, angle):
        if self.method == EULER:
            self.vel = self.vel + a(self.pos, dist, angle) * dt
            self.pos = self.pos + self.vel * dt

    def position(self,dt):
        return self.pos + self.vel* dt

    def movetime(self,dt):
        self.time = self.time + dt
        return self.time

