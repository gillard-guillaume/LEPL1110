#include <stdio.h>
#include <math.h>
#include "glfem.h"


// Guillaume Gillard
// Igor Gregoire


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    double xi[3] = {1/6, 1/6, 2/3};
    double eta[3] = {1/6, 2/3, 1/6};
    double w[3] = {1/6, 1/6, 1/6};

    
    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
        I += w[i] * f(xLoc[i], yLoc[i]);
    }




    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    


    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    double I = integrate(x,y,f);
    
//
//
//    
     
    return I;
}
