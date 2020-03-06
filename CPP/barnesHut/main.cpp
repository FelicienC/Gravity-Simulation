/* Copyright (C) Félicien Cantalloube 2020. All Rights Reserved.*/
# include "node.h"
# include "world.h"

#ifdef __APPLE__
# include <OpenGL/gl.h>
# include <OpenGL/glu.h>
# include <GLUT/glut.h>
#else
# include <GL/gl.h>
# include <GL/glu.h>
# include <GL/glut.h>
#endif


const int WIDTH = 800;
const int HEIGHT = 800;
const int NB_VERTICES = 10;
const float PI = 3.14159265;

float scale = 1e11;
float currentTime, timeBase, fps;
int frames = 0;
bool displayQuadrant = true;

World world(1000, 0);


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
  world.update();

  // Clears Color and Depth Buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  // Resets transformations
  glLoadIdentity();
  
  // Sets the camera
  gluLookAt(WIDTH/(2*scale), HEIGHT/(2*scale), 10.0f,
            WIDTH/(2*scale), HEIGHT/(2*scale),  0.0f,
            0.0f, 1.0f,  0.0f);

  // Displays all the particles
  for (int i = 0; i < world.getParticleCount(); i++) {

    Particle part = world.getParticle(i);

    // Previous trajectory
    float colorScale;
    for (int j=part.pastX.size()-1; j>=0; j--){
      colorScale = 1-j/float(part.pastX.size());
      glColor3f(part.color[0]*colorScale,
                part.color[1]*colorScale,
                part.color[2]*colorScale);
      drawCircle(part.pastX[j]/scale,
                 part.pastY[j]/scale,
                 part.radius/scale, part.nb_vertices/2);
    }

    // Current position
    glColor3f(part.color[0], part.color[1], part.color[2]);
    drawCircle(part.state[0]/scale,
               part.state[1]/scale,
               part.radius/scale, part.nb_vertices);

    // Name
    if (not part.name.empty())
      glColor3f(1, 1, 1);
      drawString((part.state[0]+part.radius)/scale,
                 (part.state[1]-part.radius)/scale,
                  0, part.name);
    }

  // Displays the quadrants
  if (displayQuadrant) {
    glColor3f(0, 1, 0);
    std::vector<Rectangle> allRects = world.getQuadrants();

    for (int j = 0; j < allRects.size(); j++) {
      glBegin(GL_LINE_STRIP);
      glVertex3f((allRects[j].x+allRects[j].size/2)/scale,
                 (allRects[j].y+allRects[j].size/2)/scale, 0.0f);
      glVertex3f((allRects[j].x-allRects[j].size/2)/scale,
                 (allRects[j].y+allRects[j].size/2)/scale, 0.0f);
      glVertex3f((allRects[j].x-allRects[j].size/2)/scale,
                 (allRects[j].y-allRects[j].size/2)/scale, 0.0f);
      glVertex3f((allRects[j].x+allRects[j].size/2)/scale,
                 (allRects[j].y-allRects[j].size/2)/scale, 0.0f);
      glVertex3f((allRects[j].x+allRects[j].size/2)/scale,
                 (allRects[j].y+allRects[j].size/2)/scale, 0.0f);
      glEnd();
    }
  }
  

  // Displays information on screen
  frames++;
  currentTime = glutGet(GLUT_ELAPSED_TIME);
  if (currentTime - timeBase > 1000) {
    fps = frames*1000.0/(currentTime-timeBase);
    timeBase = currentTime; 
    frames = 0;
  }
  glColor3f(1, 1, 1);

  drawString(-6, 4.00, 0, "Number of Objects : ");
  drawString(-4, 4.00, 0, std::to_string(world.getParticleCount()));

  drawString(-6, 3.85, 0, "Simulation Time (in days) ");
  drawString(-4, 3.85, 0, std::to_string(int(world.getSimulationTime()/86400)));

  drawString(-6, 3.70, 0, "Time step (in days) ");
  drawString(-4, 3.70, 0, std::to_string(world.getDelta()/86400));

  drawString(-6, 3.55, 0, "Frames per second : ");
  drawString(-4, 3.55, 0, std::to_string(fps));

  glutSwapBuffers();
}

/*
 * To allow zoom and escape key
 */
void processNormalKeys(unsigned char key, int x, int y) {
  
  // Escape Key
  switch (key){
    case 27:
      exit(0);
    case 'i':
      scale *= 0.9;
      break;
    case 'o':
      scale /= 0.9;
      break;
    case 'f':
      world.multiplyDelta(2.0);
      break;
    case 's':
      world.multiplyDelta(0.5);
      break;
    case 'q':
      displayQuadrant = !displayQuadrant;
  }
}


int main(int argc, char **argv) {

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
