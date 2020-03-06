#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include "node.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <valarray>
#include <deque>
#include <algorithm>


const float DELTA_MAX = 1e9;

class World
{
    public:

    World(int n, int nbtotal);

    /*
     * 
     */
    void update();

    /*
     *
     */
	void multiplyDelta(float factor);

    /*
     *
     */
	int getParticleCount() const;

    /*
     * Returns the current time in the simulation
     */
	float getSimulationTime() const;

    /*
     * returns the current time step used in the simulation
     */
	float getDelta() const;

    /*
     * Returns a specific particle
     */
	Particle getParticle(int index) const;

	/*
	 * Lists all the quadrants of the quadtree built in the simulation
	 */
	std::vector<Rectangle> getQuadrants();

    private:

    int nb_lag_total, nb_particles, nb_previous;
	float delta, simulationTime, dim;
	std::vector<Particle> particles;
	std::vector<std::valarray<double>> allVariations;
	Node Tree;

	/*
	 * Computes the squared distance between two particles
	 */
	double dist2(Particle P1, Particle P2) const;

	/*
	 * Computes the acceleration of one particle based on the distance to all 
	 * other particles. The state provided is the one of the particle wich will
	 * be accelerated.
	 */
	std::valarray<double> acceleration(std::valarray<double> oldState, int index);

	/*
	 * Implementation of the Runge-Kutta method to solve the differential 
	 * equation defined by the dynamics of the particles.
	 */
	void rungeKutta();

};

#endif // WORLD_H_INCLUDED