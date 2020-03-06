#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <valarray>
#include <deque>

const float G = 6.67e-11;
const float THETA_MAX = 1;
const float MERGE_LIMIT = 1e9;

struct Particle {
  int nb_vertices = 10;
  double mass;
  double radius;
  std::valarray<double> state = {0.0, 0.0, 0.0, 0.0};
  std::string name;
  float color[3] = {1.0f, 1.0f, 1.0f};
  std::deque<double> pastX, pastY;
};

struct Rectangle {
	float x, y, size;
};

class Node
{
  /*
   * Class representing the quadrants of the quad-tree used to increase the
   * simulation's speed.
   */
	public:

		/*
		 * Constructor using the center of the new quadrant and its size
		 */
		Node(std::valarray<double> center, float size);

		/*
		 * Recursivly enters the nodes below to compute the center of mass of
     * each of them
		 */
		std::valarray<double> computeMassCenter();

		/*
		 * Returns the force created by the node on a point, by recursively
     * dividing it until the required precision is reached.
		 */
		std::valarray<double> forceOn(std::valarray<double> point, int index);

		/*
		 * Selects the right quadrant in which to put the particle described by its
		 * x, y coordinates and size provided as input. 
		 */
		int addInCorrectQuadrant(int index, float x, float y, float r, double partMass);

	  /*
	   * Given a node, finds the quadrant below it in wich the particle fits,
     * creates a node for it and puts the particle inside
     */
		int addParticleToTree(int index, float x, float y, float r, double partMass);

		/*
		 * The memory used by the tree must be freed after each deletion of the tree
		 */
		void cleanTree();

		/*
		 * Changes the size of the root node in order to contain all the particles,
		 * event if some are porjected out
		 */
		void setSize(float newSize);

		/*
		 * Method used to draw all the quadrants of the tree
		 */
		std::vector<Rectangle> getAllQuadrants();


	private:

        Node *childNE, *childSE, *childSW, *childNW; 

        int nbParts, particleIndex;
        double mass;
        float size, particleX, particleY, particleRadius;
        std::valarray<double> center, centerOfMass; 

};

#endif // NODE_H_INCLUDED