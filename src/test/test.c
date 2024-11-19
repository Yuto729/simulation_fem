#ifndef __APPLE__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

void display () {

    /* clear window */
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-2.0,2.0,-2.0,2.0,-5.0,5.0);    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0);
    /* draw scene */
    glutWireTeapot(1);

    /* flush drawing routines to the window */
    glFlush();

}

int main ( int argc, char * argv[] ) 
{
    glutInit(&argc,argv);
    glutInitWindowSize(500,500);
    glutInitWindowPosition(0,0);
    glutInitDisplayMode(GLUT_RGB);
    glutCreateWindow("TEST OPENGL");
    glutDisplayFunc(display);
    glutMainLoop();
}