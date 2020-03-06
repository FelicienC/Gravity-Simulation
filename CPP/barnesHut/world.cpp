/* Copyright (C) Félicien Cantalloube 2020. All Rights Reserved.*/
#include "world.h"


World::World(int n, int nbtotal) :  nb_lag_total(nbtotal),
                                    nb_particles(n),
                                    nb_previous(nb_lag_total/nb_particles),
                                    delta(1e3),
                                    dim(5e11),
                                    simulationTime(0.0),
                                    Tree({0, 0}, 2*dim) {

  

  for (int i = 0; i < nb_particles; i++) {

    Particle P;
    std::valarray<double> state;
    P.mass = 1e25;
    P.radius = pow(P.mass*1e3, 0.33);
    P.state = {double(rand()%100000-50000)*1e7, // x position
               double(rand()%100000-50000)*1e7, // y position
               double(rand()%100000-50000)*1e-3, // x speed
               double(rand()%100000-50000)*1e-3}; // y speed 
    particles.push_back(P);
    allVariations.push_back(state);
  }
  // The number of points to draw is kept constant

}

void World::update(){
  int j = -1;

  std::vector<int> toDelete;
  Tree.cleanTree();
  Tree.setSize(2*dim);
  for (int i=0; i<particles.size(); i++){
     j = Tree.addParticleToTree(i,
                                particles[i].state[0],
                                particles[i].state[1],
                                particles[i].radius,
                                particles[i].mass);
     // Particle i and j should be merged together
     if (j > 0){
        // The position of the more massive particule is unchanged
        if (particles[i].mass < particles[j].mass) {
          particles[i].state[0] = particles[j].state[0];
          particles[i].state[1] = particles[j].state[1];
        }
        // Conservation of momentum
        particles[i].state[2] = particles[i].state[2] * particles[i].mass/
                                        (particles[i].mass + particles[j].mass)
                              + particles[j].state[2] * particles[j].mass/
                                        (particles[i].mass + particles[j].mass);
        particles[i].state[3] = particles[i].state[3] * particles[i].mass/
                                        (particles[i].mass + particles[j].mass)
                              + particles[j].state[3] * particles[j].mass/
                                        (particles[i].mass + particles[j].mass);
        // Mass of the new particle is the sum of both previous
        particles[i].mass = particles[i].mass + particles[j].mass;
        // Radius increases as if we had spheres not disks
        particles[i].radius = pow(particles[i].mass*1e3, 0.33);
        // The second particle is removed
        toDelete.push_back(j);
        // The variation memory is also freed
        allVariations.pop_back();
        // The total number of particle is reduced
        nb_particles--;
     }
  }
  std::sort(toDelete.begin(),toDelete.end());
  std::reverse(toDelete.begin(), toDelete.end());
  for (int i=0; i<toDelete.size(); i++){
        particles.erase(particles.begin()+toDelete[i]);
  }
  Tree.computeMassCenter();
  rungeKutta();
  simulationTime += delta;
}

int World::getParticleCount() const {
  return nb_particles;
}

float World::getSimulationTime() const {
  return simulationTime;
}

float World::getDelta() const {
  return delta;
}

Particle World::getParticle(int index) const {
  return particles[index];
}

double World::dist2(Particle P1, Particle P2) const {
  return pow(P1.state[0]-P2.state[0], 2) + pow(P1.state[1]-P2.state[1], 2);
}

void World::multiplyDelta(float factor){
  delta = (factor*delta < DELTA_MAX) ? factor*delta : delta;
}

std::valarray<double> World::acceleration(std::valarray<double> oldState,
                                          int index) {

  // dx = vx; dy = vy; dvx = ax; dvy = ay, we just need the two last ones
  std::valarray<double> state = {oldState[2], oldState[3], 0.0, 0.0};

  // We compute the acceleration as the sum of the forces ma = sum(F)
  std::valarray<double> force = Tree.forceOn({oldState[0], oldState[1]}, index);
  state[2] = force[0];
  state[3] = force[1];

  return state;
}

void World::rungeKutta() {

  for (int i = 0; i < nb_particles; i++) {

    // Computation of the 4 coefficients of the integration method
    std::valarray<double> k1 = delta * acceleration(particles[i].state, i);
    std::valarray<double> k2 = delta * acceleration(particles[i].state+k1/2, i);
    std::valarray<double> k3 = delta * acceleration(particles[i].state+k2/2, i);
    std::valarray<double> k4 = delta * acceleration(particles[i].state+k3, i);
    allVariations[i] = (k1+2*(k2+k3)+k4)/6;
  }

  // All the variations in position and speed are computed, we now apply them
  for (int i = 0; i < nb_particles; i++){

    // We add the particle position to the history
    particles[i].pastX.push_front(particles[i].state[0]);
    particles[i].pastY.push_front(particles[i].state[1]);

    // We update the position of the particle
    particles[i].state += allVariations[i];

    // We find the particle defining the dimensions of the initial quadrant
    if (abs(particles[i].state[0]) > dim)
      dim = abs(particles[i].state[0]);
    if (abs(particles[i].state[1]) > dim)
      dim = abs(particles[i].state[1]);

    // We delete the history if it contains enough positions
    if (particles[i].pastX.size() > nb_previous) {
      particles[i].pastX.pop_back();
      particles[i].pastY.pop_back();
    }
  }
}

std::vector<Rectangle> World::getQuadrants(){
  return Tree.getAllQuadrants();
}

