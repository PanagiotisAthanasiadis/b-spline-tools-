#include <iostream>
#include <fstream>
#include <vector>
#include <GL/glut.h>

using namespace std;

static long font = (long)GLUT_BITMAP_HELVETICA_18; 
static int splineOrder = 4; // Order of spline.
static float knots[12] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,90.0,100.0,110.0};
//static float knots[11] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};


// Routine to draw a bitmap character string.
void writeBitmapString(void *font,char *string)
{  
   char *c;
   for (c = string; *c != '\0'; c++) glutBitmapCharacter(font, *c);
}

float Bspline(int index, int order, float u)
{
   float coef1, coef2;
   if ( order == 1 )
   {
	  if ( index == 0 ) if ( ( knots[index] <= u ) && ( u <= knots[index+1] ) ) return 1.0;
      if ( ( knots[index] < u ) && ( u <= knots[index+1] ) ) return 1.0;
	  else return 0.0;
   }
   else
   {
      if ( knots[index + order - 1] == knots[index] ) 
	  {
	     if ( u == knots[index] ) coef1 = 1;
		 else coef1 = 0;
	  }
	  else coef1 = (u - knots[index])/(knots[index + order - 1] - knots[index]);
 
      if ( knots[index + order] == knots[index+1] )
	  {
		 if ( u == knots[index + order] ) coef2 = 1;
		 else coef2 = 0;
	  }
	  else coef2 = (knots[index + order] - u)/(knots[index + order] - knots[index+1]);
		
      return ( coef1 * Bspline(index, order-1, u) + coef2 * Bspline(index+1,order-1 ,u) );
   }
}

void specialKeyInput(int key, int x, int y)
{
   if(key == GLUT_KEY_UP) 
   {
	   if (splineOrder < 4) splineOrder++; else splineOrder = 1;
   }
   if(key == GLUT_KEY_DOWN) 
   {
	   if (splineOrder > 1) splineOrder--; else splineOrder = 4;
   }
   glutPostRedisplay();
}






void drawSpline(int index,int order)
{

    if (order == 1)
   {
	  // Spline curve.
      glBegin(GL_LINE_STRIP);
	  for (float x = knots[index]; x < knots[index+1]; x+=0.05 )
         glVertex3f( -50.0 + x, 0.0, 0.0 );	
	  glEnd();
	  glPointSize(3.0);

	  // Joints.
	  glBegin(GL_POINTS);
	  glVertex3f( -50.0 + knots[index], 0.0, 0.0);
	  glVertex3f( -50.0 + knots[index+1], 0.0, 0.0);
	  glEnd();
   }
   else {
   
    // Spline curve.
	  glBegin(GL_LINE_STRIP);
        for (float x = knots[index]; x <= knots[index + order]; x += 0.005 )
            glVertex3f( -50.0 + x, 30*Bspline(index, order, x) - 30.0, 0.0 );
            	
      glEnd();

	  // Joints.
	  glColor3f(0.0, 0.0, 0.0);
	  glBegin(GL_POINTS);
        for (int j = index; j <= index + order; j++)
            glVertex3f( -50.0 + knots[j], 30*Bspline(index, order, knots[j]) - 30.0, 0.0 );
	  glEnd();
   }
    
}


void createScene()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //clear the buffers
    glPushMatrix(); //Copy the matrix from the matrix stack ?????
    glTranslatef(0.0,0.5,0.0);



    switch (splineOrder) 
   {
      case 1:
          glRasterPos3f(-10.5, 35.0, 0.0);
          writeBitmapString((void*)font, "First-order B-splines");
	  break;
      case 2:
          glRasterPos3f(-10.5, 35.0, 0.0);
          writeBitmapString((void*)font, "Linear B-splines");
	  break;
      case 3:
          glRasterPos3f(-10.5, 35.0, 0.0);
          writeBitmapString((void*)font, "Quadratic B-splines");
	  break;
      case 4:
          glRasterPos3f(-10.5, 35.0, 0.0);
          writeBitmapString((void*)font, "Cubic B-splines");
	  break;
      default:
      break;
   }

   glEnable (GL_LINE_SMOOTH);
   glEnable (GL_BLEND);
   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
   for (int i = 0; i < 11 - splineOrder; i++ )
   {
	  switch (i) 
	  {     
         case 0:
	       glColor3f(1.0, 0.0, 0.0);
	     break;
         case 1:
	       glColor3f(0.0, 1.0, 0.0);
	     break;
         case 2:
	       glColor3f(0.0, 0.0, 1.0);
	     break;
         case 3:
	       glColor3f(1.0, 0.0, 1.0);
	     break;
         case 4:
	       glColor3f(0.0, 1.0, 1.0);
	     break;
         case 5:
	       glColor3f(1.0, 1.0, 0.0);
	     break;
         case 6:
	       glColor3f(0.0, 0.0, 0.0);
	     break;
         case 7:
	       glColor3f(1.0, 0.0, 0.0);
	     break;
         case 8:
	       glColor3f(0.5, 0.0, 0.0);
	     break;
         case 9:
	       glColor3f(0.0, 0.5, 0.0);
	     break;
         default:
         break;
	  }
       drawSpline(i, splineOrder);
   }



    glRasterPos3f(-10.5, 35.0, 0.0);
    

    //Create the x-axis
    glColor3f(0.0,0.0,0.0); //Set the drawing color to black
    glBegin(GL_LINES); //Option to draw line 
        glVertex3f(-50.0, -30.0, 0.0);
        glVertex3f(50.0, -30.0, 0.0);
    glEnd();

    //Create points for x axis
    glPointSize(6.0); //Size of the dot
    glBegin(GL_POINTS); //Point mode

        for(int i=0; i < 11; i++){glVertex3f(-50.0 + i * 10.0,-30.0,0.0);}
    glEnd();
    

    //Create the y-axis
    glBegin(GL_LINES); //Option to draw line 
        glVertex3f(-50.0, -30.0, 0.0);
        glVertex3f(-50.0, 50.0, 0.0);
    glEnd();

    //Create points for y axis
    glPointSize(6.0); //Size of the dot
    glBegin(GL_POINTS); //Point mode
        for(int i=0; i < 11; i++){glVertex3f(-50.0,-30.0 + i * 10.0,0.0);}
    glEnd();




    glPopMatrix();
    glutSwapBuffers();
}

void Render_options(int x, int u)
{
   glViewport(0, 0, (GLsizei)x, (GLsizei)u);
   glMatrixMode(GL_PROJECTION); //Defines the properties of the camera
   glLoadIdentity();//Resets us (0,0,0) for proper usage of gl0rtho
   glOrtho(-60.0, 60.0, -60.0, 60.0, -1.0, 1.0);//Define the Coordinates system
   
   glMatrixMode(GL_MODELVIEW);//Define the rules of how object are transformed(Translation,rotation,scaling)
   glLoadIdentity();//Apply the previous line of code
}

int main(int argc, char **argv)
{
    int a,b,p;//a(start),b(end),p(degree)
    vector<double> u; //knot vector
    ifstream ifile("input.txt"); //INput file data order a,b,p U{...} 
    
    if(ifile.is_open()) //File Reading 
    {
        double buf;
        ifile >> a >> b >> p ;
        while(ifile >> buf)
        {
           u.push_back(buf); 
        }
    }


    //cout << "a:" << a << " "<< "b:" << b << " " << "p:" << p << " " ;
    /*for(auto it=u.begin(); it != u.end(); it++)
    {
        cout << *it << " "; //DEBUG
    }*/

    glutInit(&argc,argv); //Initialize glut library
    glutInitDisplayMode(GLUT_DOUBLE |GLUT_RGB); // Double buffered window and rgb mode(color)
    glutInitWindowSize(800, 600); //Window size
    glutInitWindowPosition(100, 100); //Center of the screen
    glutCreateWindow("B-Spline-Demonstration"); //Name the the window
    glClearColor(1.0, 1.0, 1.0, 0.0); // set the color that clears the buffers
    glutDisplayFunc(createScene); //Call the display func in a loop
    glutReshapeFunc(Render_options); //Call the reshape func in a loop
     //glutKeyboardFunc(keyInput); //Input options
    glutSpecialFunc(specialKeyInput);
    glutMainLoop(); //Keeps the window open



    
    return 0 ;
}