#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 12:34:21 2020

@author: felicien
"""
import pygame as pg
import numpy as np

G = 6.67e-11 #m3/(kg s2)
MAX_SPEED = 0.00 #m/s
PART_MASS = 1e8  #kg
PART_RADIUS = 3 #m$
X_RANGE = 900 #m
TIME_DELTA = 100 #s

NB_PART = 20
THETA_MAX = 1

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)

class Node():
    """
    Class representing the quadrants of the quad-tree used to increase the
    simulation's speed.
    """

    def __init__(self, limit, size):
        self.mass = 0
        self.center_of_mass = None
        self.children = [None]*4 # NE, SE, SW, NW
        self.nb_parts = 0
        self.particle = None
        self.limit = limit
        self.size = size

    def compute_mass_center(self):
        """
        Recursivly enters the nodes below to compute the center of mass of
        each of them
        """

        if self.nb_parts == 1:
            self.center_of_mass = np.array([self.particle.x, self.particle.y])
            return self.center_of_mass
        self.center_of_mass = \
            sum([child.compute_mass_center()*child.mass \
                 for child in self.children if child is not None])/self.mass
        return self.center_of_mass

    def force_on(self, xy, index):
        """
        Returns the force created by the node on a point, by recursively
        dividing it until the required precision is reached.
        """

        d = (sum(map(lambda x: x**2, self.center_of_mass-xy)))**.5
        if d < 10 or (self.nb_parts == 1 and self.particle.index == index):
            return np.array([0, 0])
        if self.size/d < THETA_MAX or self.nb_parts == 1:
            return self.mass*G/d**3*(self.center_of_mass-xy)
        return sum([child.force_on(xy, index) \
                    for child in self.children if child is not None])


class Particle():
    """
    Representation of a planet
    """

    def __init__(self, index, mass, radius, x, y, vx=0, vy=0):
        self.index = index
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

    def __init__(self, center, size):

        self.particles = [Particle(i,
                                   PART_MASS,
                                   PART_RADIUS,
                                   270+np.random.randint(0, X_RANGE),
                                   np.random.randint(0, X_RANGE),
                                   2*MAX_SPEED*(np.random.rand()-.5),
                                   2*MAX_SPEED*(np.random.rand()-.5)) \
                          for i in range(NB_PART)]
        self.center = center
        self.size = size
        self.root = Node(center, size)


    def y_prime(self, t, y, index):
        """
        The time derivative of Y
        t = time at wich the derivate of y is evaluated
        y = [x, y, vx, vy]
        """

        return np.concatenate((y[2:], self.root.force_on(y[:2], index)))


    def runge_kutta(self, y_n, t_n, delta, index):
        """
        Runge Kutta integration scheme
        y_n = current value
        t_n = current time
        delta = time step
        """

        k1 = delta*self.y_prime(t_n, y_n, index)
        k2 = delta*self.y_prime(t_n+delta/2, y_n+k1/2, index)
        k3 = delta*self.y_prime(t_n+delta/2, y_n+k2/2, index)
        k4 = delta*self.y_prime(t_n+delta, y_n+k3, index)
        return y_n + (k1 + 2*(k2+k3) + k4)/6 #, t_n+delta


    def update(self, delta):
        """ Computes the new positions of the particles"""

        self.root = Node(self.center, self.size)
        for p in self.particles:
            self.add_particle(p)
        self.compute_mass()
        # Computes new positions
        for index, part in enumerate(self.particles):
            part.set_xyvxvy(self.runge_kutta(part.to_y(), 0, delta, index))


    def add_in_correct_quadrant(self, particle, current_node):
        """
        Given a node, finds the quadrant below it in wich the particle fits,
        creates a node for it and puts the particle inside
        """

        if particle.x < current_node.limit[0]:
            if particle.y < current_node.limit[1]: # SW
                if current_node.children[2] is None:
                    current_node.children[2] = \
                        Node([current_node.limit[0]-current_node.size/4,
                              current_node.limit[1]-current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[2])
            else: # NW
                if current_node.children[3] is None:
                    current_node.children[3] = \
                        Node([current_node.limit[0]-current_node.size/4,
                              current_node.limit[1]+current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[3])
        else:
            if particle.y < current_node.limit[1]: # SE
                if current_node.children[0] is None:
                    current_node.children[0] = \
                        Node([current_node.limit[0]+current_node.size/4,
                              current_node.limit[1]-current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[0])
            else: # NE
                if current_node.children[1] is None:
                    current_node.children[1] = \
                        Node([current_node.limit[0]+current_node.size/4,
                              current_node.limit[1]+current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[1])

    def add_particle(self, particle, current_node=None):
        """
        Adds a particle to the tree in the correct position, by going down from
        the root node until it reaches a empty quadrant.
        """

        # To force the start from the root
        if current_node is None:
            current_node = self.root

        # If there is nothing in that node, the particle is added there
        if current_node.nb_parts == 0:
            current_node.particle = particle
            # The particle count and the total mass is updated
            current_node.nb_parts += 1
            current_node.mass += particle.mass

        # We must subdivise the current node further to place the particle
        else:
            # The node's particle must also be brought lower if the node is a
            # leaf currently
            if current_node.nb_parts == 1:
                self.add_in_correct_quadrant(particle=current_node.particle,
                                             current_node=current_node)
                # We reset the quadrant
                current_node.particle = None
            # The particle count and the total mass of the branch is updated
            current_node.nb_parts += 1
            current_node.mass += particle.mass
            self.add_in_correct_quadrant(particle, current_node)

    def compute_mass(self):
        """
        Once all the particles are set in their quadrant, the center of mass of
        each one can be computed. This method does this.
        """

        self.root.compute_mass_center()


def main():

    # Creating the world
    world = World([270 + X_RANGE/2, X_RANGE/2], 2*X_RANGE)

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
