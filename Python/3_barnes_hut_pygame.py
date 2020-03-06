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
PART_MASS = 1e9  #kg
PART_RADIUS = 3 #m
SUN_MASS = 1e9 #kg
SUN_RADIUS = 30 #m
X_RANGE = 900 #m
TIME_DELTA = 10 #s

NB_PART = 100
THETA = 1

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)

concerned, selected = [], []

class Node():
    """
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
        if self.nb_parts == 1:
            self.center_of_mass = np.array([self.particle.x, self.particle.y])
            return self.center_of_mass
        a = sum([child.compute_mass_center()*child.mass \
                 for child in self.children if child is not None])
        self.center_of_mass = a/self.mass
        return self.center_of_mass
    
    def force_on(self, xy, index):
        d = (sum(map(lambda x:x**2, self.center_of_mass-xy)))**.5
        if d < 10:
            return np.array([0, 0])
        if self.nb_parts ==  1:
            if self.particle.index == index:
                return np.array([0, 0])
            concerned.append([self.limit[0]-self.size/2,
                              self.limit[1]-self.size/2,
                              self.size,
                              self.size])
            return self.mass*G/d**3*(self.center_of_mass-xy)
        if self.size/d < THETA:
            concerned.append([self.limit[0]-self.size/2,
                  self.limit[1]-self.size/2,
                  self.size,
                  self.size])
            return self.mass*G/d**3*(self.center_of_mass-xy)
        return sum([child.force_on(xy, index) \
                    for child in self.children if child is not None])
       

class Tree():
    """
    """
    
    def __init__(self, center, size):
        self.root = Node(center, size)
        self.all_rects = []
        self.this_rect = []
        self.centers = []
        
    def add_in_correct_quadrant(self, particle, current_node):
        if particle.x < current_node.limit[0]:
            if particle.y < current_node.limit[1]:#SW
                if current_node.children[2] is None:
                    current_node.children[2] = \
                        Node([current_node.limit[0]-current_node.size/4,
                              current_node.limit[1]-current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[2])
            else: #NW
                if current_node.children[3] is None:
                    current_node.children[3] = \
                        Node([current_node.limit[0]-current_node.size/4,
                              current_node.limit[1]+current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[3])
        else :
            if particle.y < current_node.limit[1]:#SE
                if current_node.children[0] is None:
                    current_node.children[0] = \
                        Node([current_node.limit[0]+current_node.size/4,
                              current_node.limit[1]-current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[0])
            else: #NE
                if current_node.children[1] is None:
                    current_node.children[1] = \
                        Node([current_node.limit[0]+current_node.size/4,
                              current_node.limit[1]+current_node.size/4],
                             current_node.size/2)
                self.add_particle(particle, current_node.children[1])
        
    def add_particle(self, particle, current_node=None):
        """
        """
        # To force the start from the root
        if current_node is None:
            current_node = self.root
            
        if current_node.nb_parts == 0:
            # if there is nothing in that node, the particle is added
            current_node.particle = particle
            current_node.nb_parts += 1
            current_node.mass += particle.mass
            self.all_rects.append([current_node.limit[0]-current_node.size/2,
                                   current_node.limit[1]-current_node.size/2,
                                   current_node.size,
                                   current_node.size])
        else:
            # we must subdivise the current node further to place the particle
            
            # The node's particle must also be brought lower if the node is a
            # leaf currently
            if current_node.nb_parts == 1:
                self.add_in_correct_quadrant(current_node.particle, current_node)
                current_node.particle = None
            
            current_node.nb_parts += 1
            current_node.mass += particle.mass
            self.add_in_correct_quadrant(particle, current_node)
            
    def compute_mass(self):
        self.root.compute_mass_center()
        self.centers = self.list_all_centers()
    
    def list_all_centers(self, node=None):
        if node is None:
            node = self.root
        if node.nb_parts == 1:
            return [node.center_of_mass]
        return [node.center_of_mass] +sum([self.list_all_centers(child) for child in node.children if child is not None], [])
        
    def force_on(self, xy, index):
        """ need to trace the rectangles from here"""
        return self.root.force_on(xy, index)
        
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

    def __init__(self):

        self.particles = [Particle(i,
                                   PART_MASS,
                                   PART_RADIUS,
                                   270+np.random.randint(0, X_RANGE),
                                   np.random.randint(0, X_RANGE),
                                   2*MAX_SPEED*(np.random.rand()-.5),
                                   2*MAX_SPEED*(np.random.rand()-.5)) \
                          for i in range(NB_PART)]
        
        self.quad_tree = None

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


    def y_prime(self, t, y, index):
        """
        The time derivative of Y
        t = time at wich the derivate of y is evaluated
        y = [x, y, vx, vy]
        """
        return np.concatenate((y[2:], self.quad_tree.force_on(y[:2], index)))


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
        global concerned, selected
        self.quad_tree = Tree([270 + X_RANGE/2, X_RANGE/2], 2*X_RANGE)
        for p in self.particles:
            self.quad_tree.add_particle(p)
        self.quad_tree.compute_mass()
        # Computes new positions
        for index, part in enumerate(self.particles):
            part.set_xyvxvy(self.runge_kutta(part.to_y(), 0, delta, index))
            if index == 0:
                selected = concerned[::]
            concerned = []

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
        
#        for rect in world.quad_tree.all_rects:
#            pg.draw.rect(surf, (0, 255, 0), rect, 1)
#
        for rect in selected:
            pg.draw.rect(surf, (255, 0, 0), rect, 1)
        
        pg.draw.circle(surf, (255, 0, 0), (int(world.particles[0].x), int(world.particles[0].y)), 5, 0)
            
#        for center in world.quad_tree.centers:
#            pg.draw.circle(surf, (255,0,0), (int(center[0]), int(center[1])), 1, 0)
            
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
