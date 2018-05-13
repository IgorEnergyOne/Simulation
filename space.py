#
import math
import numpy as np
import vectors as v
import unit
import planet

class Space:
    def __init__(self):
        self.planet = None
        self.unit = None


    def step(self, dt, dist):
        self.unit.move(self.acceleration, dt, dist)

    def step_rot(self, dt, dist, angle):
        self.unit.move_rot(self.acceleration, dt, dist, angle)

    def acceleration(self, pos, dist, angle):
        return self.planet.acceleration_rotational(self.unit.pos, dist,angle)


    def simulate(self, time, dist):
        max_dt = 0.1

        while(time > max_dt):
            time -= max_dt
            self.step(max_dt, dist)

        self.step(time, dist)



    
