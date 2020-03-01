#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 20:46:55 2020

@author: felicien
"""
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

G = 6.67e-11 #m3/(kg s2)
MAX_SPEED = 0.01 #m/s
PART_MASS = 1e7  #kg
PART_RADIUS = 10 #m
SUN_MASS = 1e9 #kg
SUN_RADIUS = 30 #m
X_RANGE, Y_RANGE = 1400, 900 #m
TIME_DELTA = 100 #s

NB_PART = 5


WHITE = (255, 255, 255)
BLACK = (0, 0, 0)

class Particle():
    """
    Representation of a planet
    """

    def __init__(self, mass, radius, x, y, vx=0, vy=0):
        self.mass = mass
        self.radius = radius
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

    def to_y(self):
        """ Returns the state vector of a particle [x, y, vx, vy]"""
        return [self.x, self.y, self.vx, self.vy]

    def set_xyvxvy(self, y):
        """ Sets the state of a particle according to a state vector"""
        self.x = y[0]
        self.y = y[1]
        self.vx = y[2]
        self.vy = y[3]

class World():
    """
    Simple 2D world with a central planet.
    """
    
    def __init__(self):

        self.particles = [Particle(PART_MASS,
                                   PART_RADIUS,
                                   np.random.randint(0, X_RANGE),
                                   np.random.randint(0, Y_RANGE),
                                   2*MAX_SPEED*(np.random.rand()-.5),
                                   2*MAX_SPEED*(np.random.rand()-.5)) \
                          for _ in range(NB_PART)]

        self.particles.append(Particle(SUN_MASS,
                                       SUN_RADIUS,
                                       X_RANGE/2, Y_RANGE/2, 0, 0))

    def acceleration(self, Y):
        """
        Computes the acceleration of one particle based on the distance to the
        other particles of the space
        """
        acc = np.array([0.0, 0.0])
        for part in self.particles:
            d3 = ((part.x-Y[0])**2 + (part.y-Y[1])**2)**1.5
            if d3 < 1e5: # When particles are too close, divergences occur
                continue
            acc += (part.mass/d3)*G*np.array([part.x-Y[0], part.y-Y[1]])
        return acc
    
    
    def y_prime(self, t, y):
        """ 
        The time derivative of Y
        t = time at wich the derivate of y is evaluated
        y = [x, y, vx, vy]
        """
        return np.concatenate((y[2:], self.acceleration(y)))
    
    
    def runge_kutta(self, y_n, t_n, delta):
        """
        Runge Kutta integration scheme 
        y_n = current value
        t_n = current time
        delta = time step
        """
        k1 = delta*self.y_prime(t_n, y_n)
        k2 = delta*self.y_prime(t_n+delta/2, y_n+k1/2)
        k3 = delta*self.y_prime(t_n+delta/2, y_n+k2/2)
        k4 = delta*self.y_prime(t_n+delta, y_n+k3)
        return y_n + (k1 + 2*(k2+k3) + k4)/6 #, t_n+delta
    
    def update(self, delta):
        """ Computes the new positions of the particles"""
        # Computes new positions
        for part in self.particles:
            part.set_xyvxvy(self.runge_kutta(part.to_y(), 0, delta))

# Creating the world
world = World()

trajectories = []

for i in tqdm(range(10000)):
    world.update(TIME_DELTA)
    trajectories.append([[part.x, part.y] for part in world.particles])
    
    
trajectories = np.array(trajectories)
trajectories = np.moveaxis(trajectories, 0, 1)

for particle in trajectories:
    plt.plot(particle[:,0], particle[:,1])
