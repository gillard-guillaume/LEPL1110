#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0.0;
    double xLoc[3] = {0.0, 0.0, 0.0};
    double yLoc[3] = {0.0, 0.0, 0.0};

    // magic
    double xi[3] = {1.0/6, 1.0/6, 2.0/3};
    double eta[3] = {1.0/6, 2.0/3, 1.0/6};
    double w[3] = {1.0/6, 1.0/6, 1.0/6};

    double J = fabs((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]));
    
    for (int k=0; k<3; k++){
        xLoc[k] = x[0] + xi[k]*(x[1]-x[0]) + eta[k]*(x[2]-x[0]);
        yLoc[k] = y[0] + xi[k]*(y[1]-y[0]) + eta[k]*(y[2]-y[0]);
    }

    for (int k=0; k<3; k++) I+= J*w[k]*f(xLoc[k],yLoc[k]);

  glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
  glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
  glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    return I;
}

double mitosis(double x[3], double y[3], double new_x[3], double new_y[3]){
    for (int i=0; i<3; i++){
        new_x[i] = 0.0;
        new_y[i] = 0.0;
    }
    // Matrice en Raw Major
    double P[9] = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0};
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            new_x[i] += P[3*i+j]*x[j]/2;
            new_y[i] += P[3*i+j]*y[j]/2;
        }
    }
    return 0;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n){

    if (n==1) return integrate(x,y,f);

    double new_x[3] = {0.0, 0.0, 0.0};
    double new_y[3] = {0.0, 0.0, 0.0};
    mitosis(x,y,new_x,new_y);

    double x1[3] = {x[0], new_x[0], new_x[2]};
    double y1[3] = {y[0], new_y[0], new_y[2]};

    double x2[3] = {new_x[0], x[1], new_x[1]};
    double y2[3] = {new_y[0], y[1], new_y[1]};

    double x3[3] = {new_x[2], new_x[1], x[2]};
    double y3[3] = {new_y[2], new_y[1], y[2]};

    double x4[3] = {new_x[0], new_x[1], new_x[2]};
    double y4[3] = {new_y[0], new_y[1], new_y[2]};

    double I = 0;

    I += integrateRecursive(x1,y1,f,n-1);
    I += integrateRecursive(x2,y2,f,n-1);
    I += integrateRecursive(x3,y3,f,n-1);
    I += integrateRecursive(x4,y4,f,n-1);

    return I;
    }
