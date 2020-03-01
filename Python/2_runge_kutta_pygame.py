#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 12:34:21 2020

@author: felicien
"""
import pygame as pg
import numpy as np

G = 6.67e-11 #m3/(kg s2)
MAX_SPEED = 0.01 #m/s
PART_MASS = 1e6  #kg
PART_RADIUS = 10 #m
SUN_MASS = 1e9 #kg
SUN_RADIUS = 30 #m
X_RANGE, Y_RANGE = 1400, 900 #m
TIME_DELTA = 100 #s

NB_PART = 20


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

        # Central massive particle
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
            if d3 < 1e4: # When particles are too close, divergences occur
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


def main():
    # Creating the world
    world = World()

    # Setting up pygame
    pg.init()

    # Setting up the window, font and clock
    surf = pg.display.set_mode((0, 0), pg.FULLSCREEN)
    font = pg.font.Font(None, 30)
    clock = pg.time.Clock()

    # Main game loop
    run = True
    while run:
        # Computes and updates the positions of the particles
        world.update(TIME_DELTA)

        # Draws the situation and the frame rate
        surf.fill(BLACK)
        for p in world.particles:
            pg.draw.circle(surf, WHITE, (int(p.x), int(p.y)), int(p.radius), 0)
        surf.blit(font.render(str(int(clock.get_fps())), True, WHITE), (9, 50))
        pg.display.flip()
        clock.tick()

        # Allows to quit the simulation
        for event in pg.event.get():
            if event.type == pg.QUIT:
                pg.quit()
                run = False
                break
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_ESCAPE:
                    pg.quit()
                    run = False
                    break

main()
