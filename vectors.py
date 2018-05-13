# math operations
import math
import numpy as np

def length(a): #return vector length
    return math.sqrt(a.dot(a))

def squarelength(a):
    return a.dot(a)

#def cross(a,b): #
    #c = np.array([[a[1]*b[2] - a[2]*b[1], a[0]*b[2] - a[2]*b[0], a[0]*b[1] - a[1]*b[0]])
    #return c
