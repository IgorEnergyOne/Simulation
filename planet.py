
import math
import numpy as np
import vectors as v
from constants import G

d = 0.2 # by default distance between two masses (r = 1)

class Planet:
    def __init__(self,mass = 1/G):
        self.mass = mass

    def hello(self):
        print ('Hello world!')

    def acceleration(self, pos, d):
        A = - G*(pos + np.array([d,0.0,0.0]))*(self.mass/2)/(v.length(pos + np.array([d,0.0,0.0]))**3) # first mass contibution
        B = - G*(pos + np.array([-d,0.0,0.0])) * (self.mass/2)/(v.length(pos + np.array([-d,0.0,0.0]))**3) # second mass contibution
        return A+B # full acceleration

    def acceleration_rotational(self, pos, d, angle):
        x_r = d*math.cos(angle)
        y_r = d*math.sin(angle)
        A = - G*(pos + np.array([x_r,y_r,0.0]))*(self.mass/2)/(v.length(pos + np.array([x_r,y_r,0.0]))**3) # first mass contibution
        B = - G*(pos + np.array([-x_r,-y_r,0.0])) * (self.mass/2)/(v.length(pos + np.array([-x_r,-y_r,0.0]))**3) # second mass contibution
        return A+B # full acceleration

    
