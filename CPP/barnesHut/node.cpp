/* Copyright (C) Félicien Cantalloube 2020. All Rights Reserved.*/
#include "node.h"


Node::Node(std::valarray<double> center, float squaresize) : center(center),
                                                      size(squaresize),
                                                      nbParts(0),
                                                      mass(0),
                                                      childSW(NULL),
                                                      childNW(NULL),
                                                      childNE(NULL),
                                                      childSE(NULL),
                                                      particleX(0),
                                                      particleY(0),
                                                      particleIndex(-1),
                                                      particleRadius(-1){

}

std::valarray<double> Node::computeMassCenter() {
    if(nbParts == 1) {
      centerOfMass = {particleX, particleY};
      return centerOfMass;
    }
    centerOfMass = {0.0, 0.0};
    if (childNE != NULL)
      centerOfMass += childNE->computeMassCenter()*childNE->mass;
    if (childSE != NULL)
      centerOfMass += childSE->computeMassCenter()*childSE->mass;
    if (childSW != NULL)
      centerOfMass += childSW->computeMassCenter()*childSW->mass;
    if (childNW != NULL)
      centerOfMass += childNW->computeMassCenter()*childNW->mass;
    centerOfMass /= mass;
    return centerOfMass;
}

std::valarray<double> Node::forceOn(std::valarray<double> point, int index){

    double dist = sqrt(pow(centerOfMass[0]-point[0], 2) +
                       pow(centerOfMass[1]-point[1], 2));

    if ((dist < 1e9) || (nbParts == 1 && particleIndex == index))
        return {0.0, 0.0};

    if (size/dist < THETA_MAX || nbParts == 1)
        return {mass * (centerOfMass[0]-point[0]) * G/pow(dist, 3),
                mass * (centerOfMass[1]-point[1]) * G/pow(dist, 3)};

    std::valarray<double> force = {0, 0};
    if (childNE!=NULL)
      force += childNE->forceOn(point, index);
    if (childSE!=NULL)
      force += childSE->forceOn(point, index);
    if (childSW!=NULL)
      force += childSW->forceOn(point, index);
    if (childNW!=NULL)
      force += childNW->forceOn(point, index);
    return force;
}

int  Node::addParticleToTree(int index, float x, float y, float r, double partMass) {
  
  if (nbParts == 0){
    
    // If there is nothing already in that node, the particle is added there
    particleIndex = index;
    particleX = x;
    particleY = y;
    particleRadius = r;
    // The particle count and the total mass is updated
    nbParts ++;
    mass = partMass;
    return -1;

  } else {

    // We first make sure the particle does not merge with the existing one
    if (size < 2*particleRadius+2*r){
      // We keep the position of the more massive of the two
      if(mass < partMass){
        particleX = x;
        particleY = y;
        particleRadius = r;
      }
      mass += partMass;
      return particleIndex;
    }

    // We must subdivise the current node further to place the particle

    // The node's particle must also be brought lower if the node is a
    // leaf currently (we are sure we will never get a merger here, so we don't
    // record the return value of the function)
    if (nbParts == 1){
      addInCorrectQuadrant(particleIndex, particleX, particleY, particleRadius, mass);
      particleIndex = -1;
      particleRadius = 0;
    }
  
    // we bring the wanted particle lower
    nbParts ++;
    mass += partMass;
    return addInCorrectQuadrant(index, x, y, r, partMass);
  }
    
}

int Node::addInCorrectQuadrant(int index, float x, float y, float r, double partMass) {
  std::valarray<double> newcenter;
  if (x < center[0]) {
    if (y < center[1]) { 
      if (childSW == NULL) { // South West Quadrant
        newcenter = {center[0]-size/4, center[1]-size/4};
        childSW = new Node(newcenter, size/2);
      }
      return childSW->addParticleToTree(index, x, y, r, partMass);
    } else { 
      if (childNW == NULL) { // North West Quadrant
        newcenter = {center[0]-size/4, center[1]+size/4};
        childNW = new Node(newcenter, size/2);
      }
      return childNW->addParticleToTree(index, x, y, r, partMass);
    }
  } else {
    if (y < center[1]) {
      if (childSE == NULL) { // South East Quadrant
        newcenter = {center[0]+size/4, center[1]-size/4};
        childSE = new Node(newcenter, size/2);
      }
      return childSE->addParticleToTree(index, x, y, r, partMass);
    } else {
      if (childNE == NULL) { // North East Quadrant
        newcenter = {center[0]+size/4, center[1]+size/4};
        childNE = new Node(newcenter, size/2);
      }
      return childNE->addParticleToTree(index, x, y, r, partMass);
    }
  }
}

void Node::cleanTree() {
  if (childNE != NULL){
    childNE->cleanTree();
    delete childNE;
    childNE = NULL;
  }
  if (childSE != NULL){
    childSE->cleanTree();
    delete childSE;
    childSE = NULL;
  }
  if (childSW != NULL){
    childSW->cleanTree();
    delete childSW;
    childSW = NULL;
  }
  if (childNW != NULL){
    childNW->cleanTree();
    delete childNW;
    childNW = NULL;
  }
  nbParts = 0;
  mass = 0;
}

std::vector<Rectangle> Node::getAllQuadrants() {
  // Recursively builds a list of the quadrants constructed with this tree

  std::vector<Rectangle> rectList;

  if (childNE != NULL){
    std::vector<Rectangle> a = childNE->getAllQuadrants();
    rectList.insert(rectList.begin(), a.begin(), a.end());
  }
  if (childSE != NULL){
    std::vector<Rectangle> b = childSE->getAllQuadrants();
    rectList.insert(rectList.begin(), b.begin(), b.end());
  }
  if (childSW != NULL){
    std::vector<Rectangle> c = childSW->getAllQuadrants();
    rectList.insert(rectList.begin(), c.begin(), c.end());
  }
  if (childNW != NULL){
    std::vector<Rectangle> d = childNW->getAllQuadrants();
    rectList.insert(rectList.begin(), d.begin(), d.end());
  }

  Rectangle thisRect;
  thisRect.x = center[0];
  thisRect.y = center[1];
  thisRect.size = size;
  rectList.push_back(thisRect);
  return rectList;
}

void Node::setSize(float newSize){
  size = newSize;
}