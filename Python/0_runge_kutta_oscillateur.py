#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 20:46:55 2020

@author: felicien
"""
import numpy as np
import matplotlib.pyplot as plt

def F(t, y):
    """ Damped oscillator """
    return -np.array([y[0]+0.1*abs(y[2])*y[2], y[1]+0.1*abs(y[3])*y[3]])

def y_prime(t, y):
    """
    The time derivative of Y
        t = time at wich the derivate of y is evaluated
        y = [x, y, vx, vy]
    """
    return np.concatenate((y[2:], F(t, y)))

def runge_kutta(y_prime, y_n, t_n, delta):
    """
    Runge Kutta integration scheme 
        y_prime: function (y,t)->y'
        y_n : current value
        t_n : current time
        delta : time step
    """
    k1 = delta*y_prime(t_n, y_n)
    k2 = delta*y_prime(t_n+delta/2, y_n+k1/2)
    k3 = delta*y_prime(t_n+delta/2, y_n+k2/2)
    k4 = delta*y_prime(t_n+delta, y_n+k3)
    return y_n + (k1 + 2*(k2+k3) + k4)/6, t_n+delta
    

# Initial Conditions (x, y) = (5, -3), Vx = 1
t, Y = 0, [5, -3, 1, 0]
memory = []
for i in range(500):
    Y, t = runge_kutta(y_prime, Y, t, delta=.1)
    memory.append(Y)

# Displaying the trajectory
memory = np.array(memory)
plt.plot(memory[:, 0], memory[:, 1])