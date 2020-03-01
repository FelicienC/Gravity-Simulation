# Gravity Simulation

Python and C++ simulation of gravity, using Runge-Kutta integration methods.

<p align="center">
  <img src="images/final.gif", width=400>
</p>

I started by building a small python implementation with Euler integration, but realized quickly it would not work very well (especially with pygame on top of it). So I used RK4, and switched to C++ to improve the speed.

<p align="center">
  <img src="images/pygame_21_particles.png">
</p>

It was a lot better, and looked already quite fun, even without collisions handling.
<p align="center">
  <img src="images/glut_100_particles_no_collisions.gif", width=400>
</p>

I then implemented particle fusion, as shown below.
<p align="center">
  <img src="images/glut_500_particles_merging.gif", width=400>
</p>

Added a small lag to display the previous positions, particle and fps count.
<p align="center">
  <img src="images/glut_200_particles_with.gif", width=400>
</p>

The final step was too load data from a file to look at the solar system for instance.

<p align="center">
  <img src="images/final.gif", width=400>
</p>

The next evolution would be to use a tree to sort the particles in space and turn this naive n2 approach into a Barnes-Hut simulation.
