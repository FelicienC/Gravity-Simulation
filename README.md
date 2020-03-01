# Gravity Simulation

Python and C++ simulation of gravity, using Runge-Kutta integration methods.

<p align="center">
  <img src="Images/final.gif">
</p>

I started by building a small python implementation with Euler integration, but realized quickly it would not work very well, even with small amounts of particles (especially using pygame for display):
<p align="center">
  <img src="Images/pygame_21_particles.gif">
</p>

I hence used RK4, and switched to C++ to improve the speed. It was a lot better, and looked already quite fun, even without collisions handling.
<p align="center">
  <img src="Images/glut_100_particles_no_collisions.gif">
</p>

I then implemented particle fusion, as shown below.
<p align="center">
  <img src="Images/glut_500_particles_merging.gif">
</p>

Added a small lag to display the previous positions, particle and fps count.
<p align="center">
  <img src="Images/glut_200_particles_with.gif">
</p>

The final step was too load data from a file to look at the solar system for instance.

<p align="center">
  <img src="Images/final.gif">
</p>

The next evolution would be to use a tree to sort the particles in space and turn this naive n2 approach into a Barnes-Hut simulation.
