/* Copyright (C) Félicien Cantalloube 2020. All Rights Reserved.*/

#include <stdlib.h>
#include <math.h>

#ifdef __APPLE__
# include <OpenGL/gl.h>
# include <OpenGL/glu.h>
# include <GLUT/glut.h>
#else
# include <GL/gl.h>
# include <GL/glu.h>
# include <GL/glut.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <valarray>
#include <deque>

const int WIDTH = 800;
const int HEIGHT = 800;
const int NB_VERTICES = 10;
const int MAX_POINTS = 2000;
const float INITIAL_SPEED = 0.00;
const float PI = 3.14159265;
const float G = 6.67e-11;
const float DELTA_MAX = 1e8;

struct Particle {
  double mass;
  double radius;
  std::valarray<double> state = {0.0, 0.0, 0.0, 0.0};
  std::string name;
  float color[3] = {1.0f, 1.0f, 1.0f};
  std::deque<double> pastX, pastY;
};

unsigned int seed = 42;

int nb_particles = 200;
int nb_previous;
int frames = 0;
double scale = 1e11;
double delta = 1e3;
double simulationTime = 0.0;
float currentTime, timeBase, fps;
std::vector<Particle> particles;
std::vector<std::valarray<double>> all_var;

/*
 * Computes the acceleration of one particle based on the distance to all the
 * other particles. The state provided is the one of the particle wich will
 * be accelerated.
 */
std::valarray<double> acceleration(std::valarray<double> oldState, int index) {

  // dx = vx; dy = vy; dvx = ax; dvy = ay, we just need the two last ones
  std::valarray<double> state = {oldState[2], oldState[3], 0.0, 0.0};

  // We compute the acceleration as the sum of the forces ma = sum(F)
  for (int i = 0; i < nb_particles; i++) {

    if (i != index) {

      // The force is in d square, but to project it on x and y we need the cube
      double dCube = pow(pow(particles[i].state[0]-oldState[0], 2) +
                         pow(particles[i].state[1]-oldState[1], 2), 1.5);
      // To avoid divergences with too strong forces 
      if (dCube > 1e24) {
        state[2] += particles[i].mass*G*(particles[i].state[0]-oldState[0])
                                                                         /dCube;
        state[3] += particles[i].mass*G*(particles[i].state[1]-oldState[1])
                                                                         /dCube;
      }
      
    }
    
  }

  return state;
}

/*
 * Implementation of the Runge-Kutta method to solve the differential equation
 * defined by the dynamics of the particles.
 */
void rungeKutta(double delta) {

  for (int i = 0; i < nb_particles; i++) {

    // Computation of the 4 coefficients of the method
    std::valarray<double> k1 = delta * acceleration(particles[i].state, i);
    std::valarray<double> k2 = delta * acceleration(particles[i].state+k1/2, i);
    std::valarray<double> k3 = delta * acceleration(particles[i].state+k2/2, i);
    std::valarray<double> k4 = delta * acceleration(particles[i].state+k3, i);

    // The final variation is computed with more weight for the central coefs
    all_var[i] = (k1+2*(k2+k3)+k4)/6;
  }

  // All the variations in position and speed are computed, we now apply them
  for (int i = 0; i < nb_particles; i++){
    particles[i].pastX.push_front(particles[i].state[0]);
    particles[i].pastY.push_front(particles[i].state[1]);
    particles[i].state += all_var[i];
    if (particles[i].pastX.size() > nb_previous) {
      particles[i].pastX.pop_back();
      particles[i].pastY.pop_back();
    }
  }

}

double dist2(Particle P1, Particle P2) {
  return pow(P1.state[0]-P2.state[0], 2) + pow(P1.state[1]-P2.state[1], 2);
}

/*
 * Checking for collisions between all particles in order to merge the
 * colliding ones by naively comparing n(n-1)/2 distances
 */
void mergeParticles() {

  int i = 0;
  while (i < nb_particles) {

    int j = i+1;
    while (j < nb_particles) {

      // We check for a collision
      if (dist2(particles[i], particles[j]) <
          pow(particles[i].radius + particles[i].radius, 2)) {

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
        particles.erase(particles.begin()+j);

        // The variation memory is also freed
        all_var.pop_back();

        // The total number of particle is reduced
        nb_particles--;

      } else {

        // No colision with this particle, we check the next one
        j++;
      }
    }
    i++;
  }
}

/*
 * Initializes the particle list with random positions
 */
void initParticles() {

  Particle P;
  std::valarray<double> state;
  for (int i = 0; i < nb_particles; i++) {
    P.mass = 1e25;
    P.radius = pow(P.mass*1e3, 0.33);
    P.state = {double(rand_r(&seed)%1000-500)*1e9, // x position
               double(rand_r(&seed)%1000-500)*1e9, // y position
               0.0, // x speed
               0.0}; // y speed 
    particles.push_back(P);
    all_var.push_back(state);
  }
  // The number of points to draw is kept constant
  nb_previous = MAX_POINTS/nb_particles;
}

/*
 * Reads the planets file to initialize the particles
 */
void initFromFile(char* fileName) {

  std::ifstream myfile (fileName);
  std::string line;
  std::string delimiter = ";";
  nb_particles = 0;

  if (myfile.is_open()) {
    getline (myfile,line); // Not considering the first line 
    while ( getline (myfile,line) ) {
      Particle P;
      std::valarray<double> state;
      size_t pos = 0;
      unsigned int hex_color;

      pos = line.find(delimiter);
      std::istringstream ss(line.substr(0, pos));
      ss >> P.mass;
      line.erase(0, pos + delimiter.length());

      pos = line.find(delimiter);
      std::istringstream ss2(line.substr(0, pos));
      ss2 >> P.radius;
      line.erase(0, pos + delimiter.length());

      pos = line.find(delimiter);
      std::istringstream ss3(line.substr(0, pos));
      ss3 >> P.state[0];
      line.erase(0, pos + delimiter.length());

      pos = line.find(delimiter);
      std::istringstream ss4(line.substr(0, pos));
      ss4 >> P.state[3];
      line.erase(0, pos + delimiter.length());

      pos = line.find(delimiter);
      std::istringstream ss5(line.substr(2, pos-2));
      ss5 >> std::hex >> hex_color;
      P.color[0] = ((hex_color >> 16) & 0xff) / 255.0;
      P.color[1] = ((hex_color >> 8) & 0xff) / 255.0;
      P.color[2] = (hex_color & 0xff) / 255.0;
      line.erase(0, pos + delimiter.length());

      pos = line.find(delimiter);
      std::istringstream ss6(line.substr(0, pos));
      ss6 >> P.name;
      
      nb_particles++;
      particles.push_back(P);
      all_var.push_back(state);
    }
    myfile.close();
  }
  else std::cout << "Unable to open file"; 

  nb_previous = MAX_POINTS/nb_particles;
}

void changeSize(int w, int h) { 

  // Prevent a divide by zero (you cant make a window of zero width).
  if (h <= 0)
    h = 1;
  float ratio =  w * 1.0 / h;
  // Use the Projection Matrix
  glMatrixMode(GL_PROJECTION);
  // Reset Matrix
  glLoadIdentity();
  // Set the viewport to be the entire window
  glViewport(0, 0, w, h);
  // Set the correct perspective.
  gluPerspective(45.0f, ratio, 0.1f, 100.0f);
  // Get Back to the Modelview
  glMatrixMode(GL_MODELVIEW);
}


/*
 * Draws a polygon on the screen to represent a circle
 */
void drawCircle(double x, double y, double radius, int nbPts) {

  // Prevents invisible points
  radius = (radius < 1e-2) ? 1e-2 : radius;
  nbPts = (nbPts < 3) ? 3 : nbPts;

  // Represent disk by a polygon
  glBegin(GL_POLYGON);
  for (int i = 0; i < nbPts; i++)
    glVertex3f(x+radius*cos(2*PI*i/nbPts), y+radius*sin(2*PI*i/nbPts), 0.0f);
  glEnd();
}

/*
 * Draws text on the screen
 */
void drawString(float x, float y, float z, std::string string) {
  glRasterPos3f(x, y, z);
  for (char const &c: string) {
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, c); 
    }
}


/*
 * Idle function executing the simulation and the display of the simulation
 */
void runAndRender(void) {

  // Runs the simulation
  rungeKutta(delta);
  
  // Fusion of all coliding particles
  mergeParticles();

  // Clears Color and Depth Buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  // Resets transformations
  glLoadIdentity();
  
  // Sets the camera
  gluLookAt(WIDTH/(2*scale), HEIGHT/(2*scale), 10.0f,
            WIDTH/(2*scale), HEIGHT/(2*scale),  0.0f,
            0.0f, 1.0f,  0.0f);
  
  // Displays all the particles
  for (int i = 0; i < nb_particles; i++) {

    // Previous trajectory
    for (int j=particles[i].pastX.size(); j>=0; j--){
      glColor3f(particles[i].color[0]*(1-j/float(nb_previous)),
                particles[i].color[1]*(1-j/float(nb_previous)),
                particles[i].color[2]*(1-j/float(nb_previous)));
      drawCircle((particles[i].pastX[j]-particles[i].radius)/scale,
                 (particles[i].pastY[j]-particles[i].radius)/scale,
                 particles[i].radius/scale, NB_VERTICES/2);
    }

    // Current position
    glColor3f(particles[i].color[0],
              particles[i].color[1],
              particles[i].color[2]);
    drawCircle((particles[i].state[0]-particles[i].radius)/scale,
               (particles[i].state[1]-particles[i].radius)/scale,
               particles[i].radius/scale, NB_VERTICES);

    // Name
    if (not particles[i].name.empty())
      glColor3f(1, 1, 1);
      drawString((particles[i].state[0]+particles[i].radius)/scale,
                 (particles[i].state[1]-particles[i].radius)/scale,
                  0, particles[i].name);
  }

  // Displays information on screen
  frames++;
  simulationTime += delta;
  currentTime = glutGet(GLUT_ELAPSED_TIME);
  if (currentTime - timeBase > 1000) {
    fps = frames*1000.0/(currentTime-timeBase);
    timeBase = currentTime; 
    frames = 0;
  }
  glColor3f(1, 1, 1);

  drawString(-6, 4.00, 0, "Number of Objects : ");
  drawString(-4, 4.00, 0, std::to_string(nb_particles));

  drawString(-6, 3.85, 0, "Simulation Time (in days) ");
  drawString(-4, 3.85, 0, std::to_string(int(simulationTime/86400)));

  drawString(-6, 3.70, 0, "Time step (in days) ");
  drawString(-4, 3.70, 0, std::to_string(delta/86400));

  drawString(-6, 3.55, 0, "Frames per second : ");
  drawString(-4, 3.55, 0, std::to_string(fps));

  glutSwapBuffers();
}

/*
 * To allow zoom and escape key
 */
void processNormalKeys(unsigned char key, int x, int y) {
  
  // Escape Key
  if (key == 27)
    exit(0);

  // Zoom in 
  if (key == 'i')
    scale *= 0.9;

  // Zoom out
  if (key == 'o')
    scale /= 0.9;

  // Faster 
  if (key == 'f')
    delta = (2*delta < DELTA_MAX) ? 2.0*delta : delta;

  // Slower 
  if (key == 's')
    delta *= 0.5f;
}


int main(int argc, char **argv) {


  std::string file ("-file");
  std::string nbParts ("-n");
  if (argc > 1) {
    // By default, randomly generated particles
    for (int i=0; i<argc; i++) {
      if (file.compare(argv[i]) == 0)
        initFromFile(argv[++i]);
      else if (nbParts.compare(argv[i]) == 0) {
        nb_particles = std::stoi(argv[++i]);
        initParticles();
      }
    }
   } else {
    initParticles();
   }

  // Initializes GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutCreateWindow("Gravity Simulator");

  // Switches to full screen
  glutFullScreen();

  // register callbacks
  glutDisplayFunc(runAndRender);
  glutReshapeFunc(changeSize);
  glutIdleFunc(runAndRender);
  glutKeyboardFunc(processNormalKeys);

  // enter GLUT event processing cycle
  glutMainLoop();

  return 1;
}
